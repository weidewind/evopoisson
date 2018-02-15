#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use Mutmapnolim (realdata_exists, check_realdata_restriction);
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);
use IPC::System::Simple qw(capture); 



my $protein = "h3";
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $tag = '';
my $muts = "278_INTNODE4195,209_INTNODE4241,209_INTNODE4201";
my $skip_stoppers;
my $jpeg;
GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'muts=s' => \$muts,
		'skip_stoppers' => \$skip_stoppers,
	);

my @muts = split(/,/, $muts);
my $mutmap = Mutmapnolim->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, skip_stoppers => $skip_stoppers});
my @tree_files = $mutmap->print_events(\@muts);

