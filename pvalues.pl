#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use MutMap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $config = '';
my $subtract_tallest;
my $verbose;


GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'config=s' => \$conf,
		'subtract_tallest' => \$subtract_tallest,
		'verbose'  => \$verbose,
	);

if (defined $conf && (defined $protein || defined $output || defined $input)){
	die "Mutually exclusive arguments: --config cannot be used with any other parameters";
}

## 25.01 Procedure for obtaining p-values

## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $mutmap = MutMap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, fromfile => 1});
#my @maxdepths = (50, 100, 150);
my @groups_and_names = $mutmap-> predefined_groups_and_names();
$mutmap-> concat_and_divide_simult (\@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
$mutmap-> count_pvalues(\@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}); #$self;  @restriction_levels; my @groups; my @group_names;