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
	
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, fromfile => 1}; 
print Codeversion::get_version()."\n";
sleep (180);
# hash changes
print Codeversion::get_version()."\n";
sleep (180);
print Codeversion::get_version()."\n";
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
