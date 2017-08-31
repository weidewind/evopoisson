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
use List::Util;


my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
my $no_groups;
my $verbose;
my $no_neighbour_changing;
my $no_leaves;
my $include_tips;
my $skip_stoppers;
my $mutnum_control = 0.2;
my $syn_lengths;
my $overwrite;


GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'no_groups'  => \$no_groups,
		'verbose'  => \$verbose,
		'no_neighbour_changing' => \$no_neighbour_changing,
		'no_leaves' => \$no_leaves,
		'include_tips' => \$include_tips,
		'skip_stoppers' => \$skip_stoppers,
		'mutnum_control=s' => \$mutnum_control,
		'syn_lengths' => \$syn_lengths,
		'overwrite' => \$overwrite,
	);

$| = 1;

unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest,  no_neighbour_changing => $no_neighbour_changing, no_leaves => $no_leaves, include_tips => $include_tips, skip_stoppers => $skip_stoppers,  syn_lengths => $syn_lengths, mutnum_control => $mutnum_control, fromfile => 1}; 

## Checking if appropriate realdata exists, initializing
my @restriction_levels = split(/,/, $restrictions);
my $specified_restriction = List::Util::min(@restriction_levels);
my $mutmap;
unless (Mutmapnolim::realdata_exists($args) || $overwrite) { 
	die "No realdata exists for specified parameter\n"; 
}
my $mutmap;
if ($overwrite){
	print "--overwrite option: going to overwrite existing realdata\n"; 
	$args->{fromfile} = 0;
	$mutmap = Mutmapnolim->new($args);
	$mutmap-> prepare_real_data ({restriction => $specified_restriction,step => $step});
}
else {
	my @restriction_levels = split(/,/, $restrictions);
	my $rr = Mutmapnolim::check_realdata_restriction($args);
	my $sr = List::Util::min(@restriction_levels);
	print "realdata restriction is $rr\n";
	if ($rr > $sr){ die "Error: realdata restriction is greater than minimal restriction you specified: ".$rr." > ".$sr."\n"; }
	print "Going to use existing realdata with restriction $rr\n";
	$mutmap = Mutmapnolim->new($args);
}

$mutmap-> createCodeversionFile("poisson");

my $ready = $mutmap-> count_iterations();
print "We have $ready iterations here (know nothing about their restriction, mind you)\n";

## 25.01 Procedure for obtaining p-values
my @groups_and_names;
if ($no_groups){
	@groups_and_names = $mutmap-> protein_no_group();
}
else {
	@groups_and_names = $mutmap-> predefined_groups_and_names();
}

if (!$mutnum_control){
		$mutmap-> concat_and_divide_simult_for_mutnum_controlled (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
}
else {
		$mutmap-> concat_and_divide_simult (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
}
$mutmap-> count_pvalues(\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}); #$self;  @restriction_levels; my @groups; my @group_names;

