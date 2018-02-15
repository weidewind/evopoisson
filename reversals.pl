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
my $groupfile;
my $step;
my $skip_stoppers;
GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'skip_stoppers' => \$skip_stoppers,
		'input=s' => \$input,
		'output=s' => \$output,
		'step=s' => \$step,
		'tag=s' => \$tag,
		'groupfile=s' => \$groupfile,
	);

my @group;
if ($groupfile){
	open GROUP, "<$groupfile" or die "$!\n";
	while(<GROUP>){
		my @site_node = split(/\s+/, $_);
		push @group, \@site_node;
	}
	close GROUP;
}

my $mutmap = Mutmapnolim->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, skip_stoppers => $skip_stoppers});
my %hist = $mutmap->reversals_hist($step, $tag, \@group, $groupfile);


