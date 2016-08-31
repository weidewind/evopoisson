#!/usr/bin/perl

use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use MutMap qw(realdata_exists check_realdata_restriction);
use Getopt::Long;
use Getopt::ArgvFile;
use Groups;

my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest;
my $restriction = '';
my $delete;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'restriction=i' => \$restriction,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest' => \$subtract_tallest,
		'delete'  => \$delete,
	);


## Procedure for printing real_data files
my %args = (bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest);

if  (!($overwrite) && MutMap->realdata_exists(\%args)) {
	print "Realdata with specified parameters and restriction ".check_realdata_restriction(\%args)." already exists.\n";
	print "Won't overwrite without --delete flag.\n";
	die;
}
my $mutmap = MutMap->new(\%args);
$mutmap-> prepare_real_data ($restriction);
###
