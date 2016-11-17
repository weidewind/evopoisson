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



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $subtract_tallest = '0';
my $restrictions = '50,100,150';
my $number_of_fakes = 200;
my $verbose;
my $tag = '';


GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'subtract_tallest=i' => \$subtract_tallest,
		'restrictions=s' => \$restrictions,
		'verbose'  => \$verbose,
		'number_of_fakes=i' => \$number_of_fakes,
		'tag=s' => \$tag,
	);


print "Will produce ".$number_of_fakes." fakes\n";
unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "subtract_tallest must be either 0 or 1\n";}
## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, fromfile => 1}; 
unless  (Mutmap::realdata_exists($args)) { die "No such realdata!"; }
my @restriction_levels = split(/,/, $restrictions);
my $rr = Mutmap::check_realdata_restriction($args);
my $sr = List::Util::min(@restriction_levels);
print "realdata restriction is $rr\n";
if ($rr > $sr){ die "Error: realdata restriction is greater than minimal restriction you specified: ".$rr." > ".$sr."\n"; }

## 25.01 Procedure for obtaining p-values
my $mutmap = Mutmap->new($args);
my @groups_and_names;
@groups_and_names = $mutmap-> protein_no_group();
for (my $i = 1; $i <= $number_of_fakes; $i++){
	print "Fake no $i\n";
	#my $mock_mutmap = $mutmap->mydeepclone();
	my $mock_mutmap = Mutmap->new($args);
	$mock_mutmap -> set_tag($tag."_fake_".$i); #fake realdata must be written in the subfolder
	$mock_mutmap = $mock_mutmap-> shuffle_mutator();
	#$mock_mutmap-> shuffle_mutator(); #same thing
	$mock_mutmap-> prepare_real_data ({restriction => $sr, fake => 1});
	my $newargs = {%{$args}}; # that's a reference to a new hash
	$newargs->{fake} = 1;
	$newargs->{tag} = $tag."_fake_".$i;
	$mock_mutmap = Mutmap->new($newargs); #fake realdata is taken from subfolder
	$mock_mutmap-> concat_and_divide_simult (\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
	$mock_mutmap-> count_pvalues(\@restriction_levels, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}); #$self;  @restriction_levels; my @groups; my @group_names;	
}
