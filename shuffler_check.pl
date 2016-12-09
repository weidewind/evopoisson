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



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
my $verbose;
my $tag = '';


GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'verbose'  => \$verbose,
		'tag=s' => \$tag,
	);


unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest}; 
#my @restriction_levels = split(/,/, $restrictions);

my $mutmap = Mutmap->new($args);
$mutmap->shuffle_sanity_check("_real_mutmap");
$mutmap-> prepare_real_data ({restriction => 50,step => 0.5});
$args->{fromfile} = 1;
$mutmap = Mutmap->new($args); # from file
$mutmap->shuffle_sanity_check("_real_realdata");
$mutmap = $mutmap-> shuffle_mutator();
$mutmap->shuffle_sanity_check("_shuffled");