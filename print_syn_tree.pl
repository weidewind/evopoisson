#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmapnolim;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;

my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	#

GetOptions (	
		'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		);
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest,  no_neighbour_changing => $no_neighbour_changing, no_leaves => $no_leaves, mutnum_control => $mutnum_control, fromfile => 1}; 

## Checking if appropriate realdata exists, initializing

my $mutmap;
if  (! (Mutmapnolim::realdata_exists($args))) { 
	$args->{fromfile} = 0;
	$mutmap = Mutmapnolim->new($args);
}
else {
	print "Going to use existing realdata with restriction $rr\n";
	$mutmap = Mutmapnolim->new($args);
}

$mutmap->print_tree_with_syn_lengths();