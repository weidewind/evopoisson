#!/usr/bin/perl
## launched in parallel by FDR.sh, which produces many fake samples. 
## Results are processed by grep_fdr.pl


use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;
use List::Util;
use POSIX qw(floor ceil);
use Parallel::ForkManager;
use Memory::Usage;



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
#my $number_of_fakes = 200;
my $number;
my $verbose;
my $tag = '';
my $simnumber = 10000;
my $maxmem = 4000000;
my $step = 0.5;


GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'verbose'  => \$verbose,
		#'number_of_fakes=i' => \$number_of_fakes,
		'number=i' => \$number,
		'tag=s' => \$tag,
		'simnumber=i' => \$simnumber,
		'maxmem=i' => \$maxmem,
		'step=s' => \$step,
	);


print "Will produce ".$number_of_fakes." fakes\n";
unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, fromfile => 1}; 
unless  (Mutmap::realdata_exists($args)) { die "No such realdata!"; }
my @restriction_levels = split(/,/, $restrictions);
my $rr = Mutmap::check_realdata_restriction($args);
my $sr = List::Util::min(@restriction_levels);
print "realdata restriction is $rr\n";
if ($rr > $sr){ die "Error: realdata restriction is greater than minimal restriction you specified: ".$rr." > ".$sr."\n"; }

## 25.01 Procedure for obtaining p-values

print "Fake no $number\n";
#$args->{tag} = $tag."_fake_".$number;
my $mutmap = Mutmap->new($args);
my @groups_and_names = $mutmap-> protein_no_group();
print "Fake no $number : shuffling\n";
$mutmap = $mutmap-> shuffle_mutator();
print "Fake no $number : preparing real_data\n";
 $mutmap-> prepare_real_data ({restriction => $sr, step => $step});
$args->{fake} = 1;
$args->{faketag} = $tag."_fake_".$number;
print "Fake no $number : creating mutmap from real_data\n";
$mutmap = Mutmap->new($args); #fake realdata is taken from subfolder
#$mutmap-> concat_and_divide_simult (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
#$mutmap-> count_pvalues(\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}); #$self;  @restriction_levels; my @groups; my @group_names;	

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
	unless ($max_proc_num > 0) {die "Error: Memory usage is $memusage - more than you specified with maxmem parameter\n"};
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
	      }
	    );
	$manager->run_on_finish( 
	      sub {
	         my ( $pid, $exit_code, $signal, $core ) = @_;
	         if ( $core ) {
	         	print "Process (pid: $pid) core dumped.\n";

	         } else {
	         	print "Process (pid: $pid) exited with code $exit_code and signal $signal.\n";

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
	      $manager->start and next;
	      system( $command );
	      $manager->finish;
	   }
	$manager->wait_all_children;
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

$mutmap-> concat_and_divide_simult (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
$mutmap-> count_pvalues({restriction_levels => \@restriction_levels, groups => \@{$groups_and_names[0]}, group_names => \@{$groups_and_names[1]}}); #$self;  @restriction_levels; my @groups; my @group_names;


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
	my $command = $exports.$perlocation." iterations_gulp.pl --protein $protein --state $state --subtract_tallest $subtract_tallest --iterations $its --tag $tag --restriction $sr ";
	$command = $command."--faketag ".$args->{faketag}." ";
	if($output){$command = $command."--output $output ";}
	if($input){$command = $command."--input $input ";}
	if ($verbose){ $command = $command." --verbose ";}
	if ($memusage){ $command = $command." --memusage ";}
	if ($no_leaves){ $command = $command." --no_leaves ";}
	if ($no_neighbour_changing){ $command = $command." --no_neighbour_changing ";}
	return $command;
	
}

sub mockcomm {
	my $tag = shift;
	my $command = "sleep 10";
	return $command;
}


#delete files
my $path = Mutmap::pathFinder($args);
my $unreadable = Mutmap::temp_tag();
my $victim = File::Spec->catdir($path, $unreadable);
opendir (TEMP, $victim);
my @files = readdir(TEMP);
closedir(TEMP);
foreach my $f(@files){
	my $filepath = File::Spec->catfile($victim, $f);
	unlink $filepath or warn "Cannot unlink $filepath: $!\n";
}

