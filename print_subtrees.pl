#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use Mutmap (realdata_exists, check_realdata_restriction);
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);



my $protein = "h3";
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $muts = "278,INTNODE4195,209,INTNODE4241,209,INTNODE4201";
GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'muts=s' => \$muts,
	);

my @muts = split(/,/, $muts);
for (my $i = 0; $i < scalar @muts; $i++){
	unless ($muts[$i] =~ /^[1-9]\d*$/) {die "Please, enter comma-separated sites and corresponding nodes. Example: 278,INTNODE4195,209,INTNODE4241,209,INTNODE4201";}
	$i++; #need to check only odd elements
}
my $mutmap = Mutmap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state});

$mutmap->print_subtree_with_mutations(\@muts);

