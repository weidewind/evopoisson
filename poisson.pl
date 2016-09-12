#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use MutMap;
use Locker;
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
my $maxmem = 4000000;


my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restriction = '50,100,150';
my $simnumber = 10000;
my $maxprocs = 2;
my $verbose;


GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'simnumber=i' => \$simnumber,
		'maxprocs=i' => \$maxprocs,
		'verbose'  => \$verbose,
	);



unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, fromfile => 1}; 

## Checking if appropriate realdata exists
my @restriction_levels = split(/,/, $restrictions);
my $specified_restriction = List::Util::min(@restriction_levels);

if  (! (MutMap::realdata_exists($args))) { 
	print "No realdata exists for specified parameters, going to prepare it.."; 
	$args->{fromfile} = 0;
	my $mutmap = MutMap->new($args);
	$mutmap-> prepare_real_data ($specified_restriction);
}
else {
	my $rr = MutMap::check_realdata_restriction($args);
	if ($rr > $specified_restriction){
		print "Existing realdata restriction is greater than minimal restriction you specified: ".$rr." > ".$specified_restriction."\nGoing to overwrite realdata..\n"; 
		$args->{fromfile} = 0;
		my $mutmap = MutMap->new($args);
		$mutmap-> prepare_real_data ($specified_restriction);
	}
	else {
		print "Going to use existing realdata with restriction $rr\n";
	}
}

$mu->record('just before mutmap creation');
##
## Computing the number of iterations and size of gulps
my $mutmap = MutMap->new($args); # from file
$mu->record('mutmap created');
my $ready = $mutmap-> count_iterations();
print "Already have $ready iterations (know nothing about their restriction, mind you)\n";
my $newtag = $mutmap-> iterations_maxtag() + 1;
print "New iteration tags will start from $newtag\n";

my $sim = $simnumber-$ready;
if ($sim > 0){
	my $its_for_proc = List::Util::min(500, int($sim/$maxprocs));
	if ($its_for_proc == 0) {$its_for_proc = 1};
	my $proc_num = int($sim/$its_for_proc);
	print "At least $proc_num files, each of them containing $its_for_proc iterations, will be produced\n";
	##
	## Preparing commands for Forkmanager
	my @commands;
	for (my $tag = $newtag; $tag < $proc_num+$newtag; $tag++){
		my $command = mockcomm($tag, $its_for_proc);
		push @commands, $command;
		print $command."\n";	
	}
	my $remainder = $sim%$its_for_proc;
	print "Remainder is $remainder\n";
	if ( $remainder > 0) { 
		my $command = mockcomm($proc_num+$newtag, $remainder);
		push @commands,$command;
		print $command."\n";
	}
	##
	## Forkmanager setup
	my $manager = new Parallel::ForkManager($maxprocs);
	my $locker = Locker->new($mutmap, $maxmem);

	$manager->run_on_finish( 
	      sub {
	         my ( $pid, $exit_code, $signal, $core ) = @_;
	         print "entered finish\n";
	         if ( $core ) {
	         	$locker->delete_entry();
	         	print "Process (pid: $pid) core dumped.\n";
	            $mu->record("Process (pid: $pid) core dumped.\n");
	         } else {
	         	print "Going to finish $pid..\n";
				$locker->delete_entry();
	            $mu->record("Process (pid: $pid) exited with code $exit_code and signal $signal.\n");
	         }
	      }
	   );
	$manager->run_on_wait( 
	      sub {
	         print "Waiting ... \n";
	      },
	      5 # time interval between checks
	   ); 
	##
	## Launching a series of new iteration_gulps if needed 
	my $pid;
	foreach my $command (@commands) {
		
		my $success;
		while(!$success){
		    	$success = 1;
		    	print "my pid is $pid\n"; # previous child pid
		    	sleep(10);
		    	print "I am awake\n";
		    	sleep(10);
		    	print "I am awake\n";
		    	$manager->wait_children();
		    	sleep(10);
		    	print "I am awake\n";
		    	sleep(10);
		    	print "Getting up\n";
				#$success = $locker->try_and_start(); # commented out for debugging
		 }
	  #  print "Starting $command\n";
      #  $mu->record("Starting $command\n");
	    #$pid = $manager->start && next; # fork! parent gets child pid, child gets 0. && means that "next" is never read if first argument evaluates to 0
	    $pid = $manager->start and next;
	    system( $command );
	    print (" Went out\n"); # must be 0
		$manager->finish;
	};
	$manager->wait_all_children;
	$mu->dump();
	##
}




## 25.01 Procedure for obtaining p-values
#my $mutmap = MutMap->new($args);
#my @groups_and_names = $mutmap-> predefined_groups_and_names();
#$mutmap-> concat_and_divide_simult (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
#$mutmap-> count_pvalues(\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}); #$self;  @restriction_levels; my @groups; my @group_names;


sub mycomm {
	my $tag = shift;
	my $its = shift;
	my $command = "perl iterations_gulp.pl --protein $protein --state $state --subtract_tallest $subtract_tallest --iterations $its --tag $tag ";
	if($output){$command = $command."--output $output ";}
	if($input){$command = $command."--input $input ";}
	if ($verbose){ $command = $command." --verbose ";}
	return $command;
	
}

sub mockcomm {
	my $tag = shift;
	sleep(10);
	print "$tag is OK\n";
	return "echo Hi!\n";
}