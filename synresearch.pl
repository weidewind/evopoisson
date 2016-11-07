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
		use Data::Dumper;



my $protein;
my $state = 'syn';
my $input = '';
my $output = '';	# option variable with default value
my $verbose;


GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'verbose'  => \$verbose,
	);

	my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state}; 
	my $mutmap = Mutmap->new($args);
#	my $hash = $mutmap->synmut_types();
	$mutmap->reversals_list();
	print Dumper $hash;