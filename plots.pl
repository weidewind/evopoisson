#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use Mutmapnolim;
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $verbose;
my $step;
my $skip_stoppers;
my $cumulative;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'step=s' => \$step,
		'verbose'  => \$verbose,
		'skip_stoppers' => \$skip_stoppers,
		'cumulative' => \$cumulative,
	);

my $mutmap = Mutmapnolim->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, skip_stoppers => $skip_stoppers});

## prints files for drawing the chosen version of plots: mutations_density_in_circles(distance from ancestral mutation), points correspond to mutations
## There could be more than one mutation at the given distance, so number of points on the plot can be less than the number of mutations in the subtree under analysis 
$mutmap->egor_diff_rings_site_entrenchment($step, $cumulative);





