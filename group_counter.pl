#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;
use List::Util;



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
my $verbose;


GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'verbose'  => \$verbose,
	);



unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, fromfile => 1}; 
unless  (Mutmap::realdata_exists($args)) { die "No such realdata!"; }
my @restriction_levels = split(/,/, $restrictions);
my $rr = Mutmap::check_realdata_restriction($args);
my $sr = List::Util::min(@restriction_levels);
print "realdata restriction is $rr\n";
if ($rr > $sr){ die "Error: realdata restriction is greater than minimal restriction you specified: ".$rr." > ".$sr."\n"; }

## 25.01 Procedure for obtaining counts
my $mutmap = Mutmap->new($args);
my @groups_and_names = $mutmap-> predefined_groups_and_names();
$mutmap-> group_counter(\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
