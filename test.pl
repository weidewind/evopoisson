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
use Storable qw(store retrieve lock_retrieve);
use Codeversion;


my $protein = 'h1';
my $state = 'nsyn';
my $input = '';
my $output = 'debugging_halfbin';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
my $verbose;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'verbose'  => \$verbose,
	);
	
#my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, fromfile => 1}; 
#print Codeversion::get_version()."\n";
#sleep (180);
# hash changes
#print Codeversion::get_version()."\n";

# all sites from Groups_and_sites restriction 50, sorted by adapt_median
my @h1 = (171,184,202,287,4,90,73,151,88,179,201,415,238,269,232,418,326,173,488,289,361,3,470,13,277,186,111,145,182,156,171,489,155,106,60,102,15,86,361,261,91,240,145,202,209,202,178,293,128,292,99,138,184,138,97,111,16,176,240,210,206,311,205,238,205,74,210,496,11,78,113,232,71,85,142,88,91,155,470,209,169,97,223,207,169,12,206,289,13,144);
my @h3 = (243,213,206,70,176,138,276,212,179,149,161,18,153,264,94,15,11,188,466,160,162,315,546,41,66,278,142,175,209,147,137,158,495,47,469,16,205,172,468,292,294,140,147,99,172,140,294,73,174,402,99,291,258,241,19,14,202,171,188,188,110,66,204,233,213,205,99,160,208,69,98,363,78,149,79,160,557,171,175,242,161,209,291,91,159,174,218,78,229,205,153,391,151,173,189);
my @n1 = (254,369,455,93,64,84,430,15,105,14,52,332,70,173,434,386,232,390,388,263,16,427,396,336,369,388,23,34,248,105,254,454,220,101,57,352,366,248,222,250,270,248,367,274,163,263,86,382,432,270,339,200,332,435,80,287,250,388,95,344,307,34);
my @n2 = (172,141,368,308,18,47,369,358,370,290,46,52,437,153,381,82,43,248,127,265,313,20,367,344,434,143,126,253,463,401,331,93,431,400,56,307,69,208,220,149,302,267,328,334,339,356,336,390,307,197,143,221,369,194,346,43,347,249,329,26,197,466,93,248,155,199,344,370,435,336,431,385);
foreach my $site(@h3){
	Groups::describe_site("h3", $site);
	print "\n";
}

#my $output_base = Mutmap::pathFinder ($args);
#my $realdatapath = File::Spec->catfile($output_base, $args->{protein}."_". Mutmap::state_tag($args->{state})."_realdata");
#my $realdata = lock_retrieve ($realdatapath) or die "Cannot retrieve ".$realdatapath;	
#print " step ".$realdata->{step};
#my $mutmap = Mutmap->new($args);

#my @groups_and_names = $mutmap-> predefined_groups_and_names();
#for(my $i = 0; $i <  scalar @{$groups_and_names[0]}; $i++){
#	print $groups_and_names[1][$i]."\n";
#	foreach my $site (@{$groups_and_names[0][$i]}){
#		print $site."\t";
#	}
#	print "\n";
#	print "-------\n";
#}
