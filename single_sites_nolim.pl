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



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
my $verbose;
my $tag = '';
my $mutnum_control = 0;
my $skip_stoppers;
my $no_neighbour_changing;
my $overwrite;



GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'verbose'  => \$verbose,
		'tag=s' => \$tag,
		'mutnum_control=s' => \$mutnum_control,
		'no_neighbour_changing' => \$no_neighbour_changing,
		'skip_stoppers' => \$skip_stoppers,
		'overwrite' => \$overwrite,
	);



unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, no_neighbour_changing => $no_neighbour_changing, skip_stoppers => $skip_stoppers, mutnum_control => $mutnum_control, fromfile => 1}; 
unless  (Mutmapnolim::realdata_exists($args) || $overwrite) { die "No such realdata!"; }
my $mutmap;
my @restriction_levels = split(/,/, $restrictions);
my $sr = List::Util::min(@restriction_levels);
if ($overwrite){
	print "--overwrite option: going to overwrite existing realdata\n"; 
	$args->{fromfile} = 0;
	$mutmap = Mutmapnolim->new($args);
	$mutmap-> prepare_real_data ({restriction => $sr,step => $step});
}
else {
	my $rr = Mutmapnolim::check_realdata_restriction($args);
	print "realdata restriction is $rr\n";
	if ($rr > $sr){ die "Error: realdata restriction is greater than minimal restriction you specified: ".$rr." > ".$sr."\n"; }
	$mutmap = Mutmapnolim->new($args);
}

## 25.01 Procedure for obtaining p-values for every ancestor node (site_node)
$mutmap -> set_tag($tag);
$mutmap-> createCodeversionFile("single_sites");
$mutmap-> concat_and_divide_simult_single_sites (\@restriction_levels, $mutnum_control);
$mutmap-> count_single_site_pvalues({restriction_levels => \@restriction_levels, stattypes => ["mean", "median", "bp"]}); #$self;  @restriction_levels; my @groups; my @group_names;

