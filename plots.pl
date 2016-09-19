#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use Mutmap;
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $verbose;
GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'verbose'  => \$verbose,
	);

my $mutmap = Mutmap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state});

## prints files for drawing the chosen version of plots: mutations_density_in_circles(distance from ancestral mutation), points correspond to mutations
## There could be more than one mutation at the given distance, so number of points on the plot can be less than the number of mutations in the subtree under analysis 
$mutmap->egor_smart_site_entrenchment($verbose);




