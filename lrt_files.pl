#!/usr/bin/perl
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use MutMap;
use Getopt::Long;
use Getopt::ArgvFile;



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,

	);

	
# Procedure for printing _for_LRT files
# Reads data/_distance_matrix.csv prepared by Distances.R from RCoevolution 

my $mutmap = MutMap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state});
$mutmap->set_distance_matrix();
$mutmap->print_data_for_LRT();
	