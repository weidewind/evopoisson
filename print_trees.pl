#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use MutMap;
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $sites;
GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'sites=s' => \$sites,
	);

my @sites = split(/,/, $sites);
foreach my $site(@sites){
	unless ($site =~ /^[1-9]\d*$/) {die "Please, enter comma-separated sites. Example: 23,124,98";}
}
my $mutmap = MutMap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state});
foreach my $site(@sites){
	$mutmap->print_static_tree_with_mutations($site);
}
