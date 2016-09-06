#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use MutMap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;
use List::Util;
use POSIX qw(floor ceil);
use Parallel::ForkManager;

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
##
## Counting existing number of iterations, launching a series of new iteration_gulps if needed
my $mutmap = MutMap->new($args); # from file
my $ready = $mutmap-> count_iterations();
print "Already have $ready iterations (know nothing about their restriction, mind you)\n";
my $newtag = $mutmap-> iterations_maxtag() + 1;
print "New iteration tags will start from $newtag\n";

my $sim = $simnumber-$ready;
my $its_for_proc = List::Util::min(500, int($sim/$maxprocs));
my $proc_num = int($sim/$its_for_proc);
my @commands;
for (my $tag = $newtag; $tag < $proc_num+$newtag; $tag++){
	my $command = mycomm($tag, $its_for_proc);
	push @commands, $command;
	print $command."\n";	
}
my $remainder = $sim%$its_for_proc;
print "Remainder is $remainder\n";
if ( $remainder > 0) { 
	my $command = mycomm($proc_num+$newtag, $remainder);
	push @commands,$command;
	print $command."\n";
}
#my $manager = new Parallel::ForkManager($maxprocs);
#foreach my $command (@commands) {
#      $manager->start and next;
#      system( $command );
#      $manager->finish;
#   };


## 25.01 Procedure for obtaining p-values
#my $mutmap = MutMap->new($args);
#my @groups_and_names = $mutmap-> predefined_groups_and_names();
#$mutmap-> concat_and_divide_simult (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
#$mutmap-> count_pvalues(\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}); #$self;  @restriction_levels; my @groups; my @group_names;

sub mycomm {
	my $tag = shift;
	my $its = shift;
	my $command = "perl iterations_gulp.pl -p $protein -o $output --input $input -s $state --iterations $its --tag $tag --subtract_tallest $subtract_tallest";
	if ($verbose){ $command = $command." --verbose";}
	return $command;
	
}