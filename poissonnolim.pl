#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmapnolim;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;
use List::Util;
use POSIX qw(floor ceil);
use Parallel::ForkManager;
use Memory::Usage;

my $mu = Memory::Usage->new();
$mu->record('starting work');

#

my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
my $step = 0.5;
my $simnumber = 10000;
my $maxmem = 4000000;
my $no_groups;
my $verbose;
my $no_neighbour_changing;
my $no_leaves;
my $include_tips;
my $skip_stoppers;
my $switch;
my $mutnum_control = 0.2;
my $fake;

my $lifetime_restr;
my $onestrip;
my $shuffler_type = "exp";
my $debugmode;
my $poisson; # exp shuffler is poisson not exp :)
my $syn_lengths;


GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'simnumber=i' => \$simnumber,
		'maxmem=i' => \$maxmem,
		'step=s' => \$step,
		'no_groups'  => \$no_groups,
		'verbose'  => \$verbose,
		'no_neighbour_changing' => \$no_neighbour_changing,
		'no_leaves' => \$no_leaves,
		'include_tips' => \$include_tips,
		'skip_stoppers' => \$skip_stoppers,
		'switch' =>\$switch,
		'mutnum_control=s' => \$mutnum_control,
		'fake' => \$fake,
		'lifetime_restr' => \$lifetime_restr,
		'onestrip' => \$onestrip,
		'shuffler_type=s' => \$shuffler_type,
		'debugmode' => \$debugmode,
		'distrpoisson' =>\$poisson,
		'syn_lengths' => \$syn_lengths,
	);

$| = 1;

unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest,  no_neighbour_changing => $no_neighbour_changing, no_leaves => $no_leaves, include_tips => $include_tips, skip_stoppers => $skip_stoppers, distrpoisson => $poisson, syn_lengths => $syn_lengths, mutnum_control => $mutnum_control, fromfile => 1}; 

## Checking if appropriate realdata exists, initializing
my @restriction_levels = split(/,/, $restrictions);
my $specified_restriction = List::Util::min(@restriction_levels);
my $mutmap;
if  (! (Mutmapnolim::realdata_exists($args))) { 
	print "No realdata exists for specified parameters, going to prepare it.."; 
	$args->{fromfile} = 0;
	$mutmap = Mutmapnolim->new($args);
	$mutmap-> prepare_real_data ({restriction => $specified_restriction,step => $step});
}
else {
	my $rr = Mutmapnolim::check_realdata_restriction($args);
	if ($rr > $specified_restriction){
		print "Existing realdata restriction is greater than the minimal restriction you specified: ".$rr." > ".$specified_restriction."\nGoing to overwrite realdata..\n"; 
		$args->{fromfile} = 0;
		$mutmap = Mutmapnolim->new($args);
		$mutmap-> prepare_real_data ({restriction => $specified_restriction,step => $step});
	}
	else {
		print "Going to use existing realdata with restriction $rr\n";
		$mutmap = Mutmapnolim->new($args);
	}
}

if ($fake){
	$mutmap = $mutmap-> shuffle_mutator();
	print "Fake : preparing real_data\n";
	$mutmap-> prepare_real_data ({restriction => $specified_restriction, step => $step});
}


$mu->record('just before mutmap creation');
$args->{fromfile} = 1;
my $mutmap = Mutmapnolim->new($args); # from file
$mu->record('mutmap created');
$mutmap-> createCodeversionFile("poisson");
##
## 
## Computing gulp sizes and number of concurrent processes
my $ready = $mutmap-> count_iterations();
print "Already have $ready iterations (know nothing about their restriction, mind you)\n";
my $newtag = $mutmap-> iterations_maxtag() + 1;
my $sim = $simnumber-$ready;
if ($sim > 0){
	print "New iteration tags will start from $newtag\n";
	## First iterations_gulp runs separately - to estimate memusage
	my $probe = 5; # first iterations-Gulp size, used for memusage esttimation, is not used afterwards
	my $command = mycomm("probe", $probe, "memusage"); # prints memusage in file
	system( $command );
	my $locker = Memusage->get_locker($mutmap);
	my $memusage = $locker->get_memusage();
	print "Memusage is ".$memusage."\n";
	##
	my $max_proc_num = int($maxmem/$memusage);
	print " Max proc num $max_proc_num\n";
	unless ($max_proc_num > 0) {die "Error: Memory usage is $memusage - more than you specified with masmem parameter\n"};
	my @iters;
	my $its_for_proc;
	if ($sim/$max_proc_num < 1) {
		$its_for_proc = 1;
		for (my $i = 0; $i < $sim; $i++){
			push @iters, $its_for_proc;
		}
	}
	else {
		$its_for_proc = int($sim/$max_proc_num);
		print " its for proc $its_for_proc\n";
		my $remainder = $sim%$max_proc_num;
		print " remainder $remainder\n";
		for (my $i = 0; $i < $max_proc_num; $i++){
			push @iters, $its_for_proc;
		}
		for (my $i = 0; $i < $remainder; $i++){
			$iters[$i] = $iters[$i] + 1;
		}
	}
	my $numfiles = scalar @iters;
	print "$numfiles files, each of them containing $its_for_proc or ".($its_for_proc+1)." iterations, will be produced\n";
	##
	## Preparing commands for Forkmanager
	my @commands;
	my $tag = $newtag;
	foreach my $it(@iters){	
		my $command = mycomm($tag, $it);
		push @commands, $command;
		print $command."\n";	
		$tag++;
	}
	##
	## Forkmanager setup
	my $manager = new Parallel::ForkManager(30);
	#my $lockfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}, "lock");

	$manager->run_on_start( 
	      sub {
	      	my $pid = shift;
	      	print "Starting child processes under process id $pid\n";
	        $mu->record("Process (pid: $pid) started.\n");
	      }
	    );
	$manager->run_on_finish( 
	      sub {
	         my ( $pid, $exit_code, $signal, $core ) = @_;
	         if ( $core ) {
	         	print "Process (pid: $pid) core dumped.\n";
	            $mu->record("Process (pid: $pid) core dumped.\n");
	         } else {
	         	print "Process (pid: $pid) exited with code $exit_code and signal $signal.\n";
	            $mu->record("Process (pid: $pid) exited with code $exit_code and signal $signal.\n");
	         }
	      }
	   );
	$manager->run_on_wait( 
	      sub {
	         print "Waiting for all children to terminate... \n";
	      },
	      180 # time interval between checks
	   ); 
	##
	## Launching a series of new iteration_gulps if needed 

	foreach my $command (@commands) {
		  sleep(3); # because the new server somehow managed to launch 5 procs with the same seed ) 
	      $manager->start and next;
	      system( $command );
	      $manager->finish;
	   }
	$manager->wait_all_children;
	$mu->dump();
	##
}
## 25.01 Procedure for obtaining p-values
my @groups_and_names;
if ($no_groups){
	@groups_and_names = $mutmap-> protein_no_group();
}
else {
	@groups_and_names = $mutmap-> predefined_groups_and_names();
}

$mutmap-> concat_and_divide_simult_for_mutnum_controlled (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
$mutmap-> count_pvalues(\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}); #$self;  @restriction_levels; my @groups; my @group_names;


sub mycomm {
	my $tag = shift;
	my $its = shift;
	my $memusage = shift;
	my $perlocation = "perl";
	my $exports = "";
	if ($switch) {
		$perlocation = "~/perl5/perlbrew/perls/perl-5.22.1/bin/perl";
	 	$exports = "export PERL5LIB=/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/.perl/lib/perl5/5.22.1/x86_64-linux:/export/home/popova/.perl/lib/perl5/5.22.1:/export/home/popova/.perl/lib/perl5/x86_64-linux:/export/home/popova/.perl/lib/perl5:/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/perl5/lib/perl5/x86_64-linux:/export/home/popova/perl5/lib/perl5:/export/home/popova/workspace/evopoisson:$PERL5LIB; ";
	}
	my $command = $exports.$perlocation." iterations_gulp.pl --protein $protein --state $state --subtract_tallest $subtract_tallest --iterations $its --tag $tag --restriction $specified_restriction --shuffler_type $shuffler_type ";
	if($output){$command = $command."--output $output ";}
	if($input){$command = $command."--input $input ";}
	if ($verbose){ $command = $command." --verbose ";}
	if ($memusage){ $command = $command." --memusage ";}
	if ($no_leaves){ $command = $command." --no_leaves ";}
	if ($include_tips) {$command = $command." --include_tips ";}
	if ($skip_stoppers) {$command = $command." --skip_stoppers ";}
	if ($no_neighbour_changing){ $command = $command." --no_neighbour_changing ";}
	if ($lifetime_restr) { $command = $command." --lifetime_restr ";}
	if ($onestrip) { $command = $command." --onestrip ";}
	if ($debugmode) { $command = $command." --debugmode ";}
	if ($poisson) { $command = $command." --distrpoisson ";}
	if ($shuffler_type ne "strip" && $onestrip){print "onestrip option is ignored (only meaningful for strip shuffler_type)";}
	return $command;
	
}

sub mockcomm {
	my $tag = shift;
	my $command = "sleep 10";
	return $command;
}