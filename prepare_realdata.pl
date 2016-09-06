#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use MutMap;
use Getopt::Long;
use Getopt::ArgvFile;
use Groups;

my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = 0;
my $restriction = 50;
my $delete;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'restriction=i' => \$restriction,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'delete'  => \$delete,
	);



unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1, got $subtract_tallest \n";}
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest};

if  (!($overwrite) && MutMap::realdata_exists($args)) {
	print "Checking existing realdata restriction..\n";
	print "Realdata with specified parameters and restriction ".MutMap::check_realdata_restriction($args)." already exists.\n";
	print "Won't overwrite without --delete flag.\n";
	die;
}
## Procedure for printing real_data files
my $mutmap = MutMap->new($args);
$mutmap-> prepare_real_data ($restriction);
###
