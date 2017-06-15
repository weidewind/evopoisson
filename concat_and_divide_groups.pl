#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmapnolim;
use Getopt::Long;
use Getopt::ArgvFile;
use Groups;



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $groupnumber;
my $subtract_tallest = '0';
my $verbose;
my $no_neighbour_changing;
my $no_leaves;
my $restrictions = '50,100,150';
my $faketag;

my $lifetime_restr;
my $onestrip;
my $debugmode;
my $syn_lengths;
my $poisson;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'groupnumber=i' => \$groupnumber,
		'subtract_tallest=i' => \$subtract_tallest,
		'verbose'  => \$verbose,
		'no_neighbour_changing' => \$no_neighbour_changing,
		'no_leaves' => \$no_leaves,
		'restrictions=s' => \$restrictions,
		'faketag=s' => \$faketag,
		'lifetime_restr' => \$lifetime_restr,
		'onestrip' => \$onestrip,
		'debugmode' => \$debugmode,
		'distrpoisson' =>\$poisson,
		'syn_lengths' =>\$syn_lengths,
	);


## Procedure for launching a gulp of iterations
unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "--subtract_tallest must be either 0 or 1\n";}
unless (defined $tag){die "There will be several gulp files in one folder, therefore a --tag for each gulp must be specified (and it should be an integer)";}
unless (defined $groupnumber){die "NO group number!\n";}
my @restriction_levels = split(/,/, $restrictions);
my $specified_restriction = List::Util::min(@restriction_levels);
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, no_neighbour_changing => $no_neighbour_changing, no_leaves => $no_leaves, syn_lengths => $syn_lengths, fromfile => 1, faketag => $faketag}; 
unless  (Mutmapnolim::realdata_exists($args)) { die "No such realdata!"; }
$mutmap = Mutmapnolim->new($args);
my $rr = $mutmap->{realdata}{restriction};
unless (defined $rr) {die "realdata restriction is not defined!\n";}
unless ($rr <= $specified_restriction){
		die  "Existing realdata restriction is greater than the minimal restriction you specified: ".$rr." > ".$specified_restriction."\nGoing to overwrite realdata..\n"; 
}
print "Using existing realdata with restriction $rr\n";

my @groups_and_names;
if ($no_groups){
	@groups_and_names = $mutmap-> protein_no_group();
}
else {
	@groups_and_names = $mutmap-> predefined_groups_and_names();
}
my @group_and_name;
push @{$group_and_name[0]}, $groups_and_names[0][$groupnumber];
push @{$group_and_name[1]}, $groups_and_names[1][$groupnumber];

if ($verbose) { print "Starting concat_and_divide for group $groupnumber\n"; }
$mutmap-> concat_and_divide_simult_for_mutnum_controlled (\@restriction_levels, \@{$group_and_name[0]}, \@{$group_and_name[1]});
if ($verbose) { print "Finished concat_and_divide for group $groupnumber".$groups_and_names[1][$groupnumber]."\n"; }
###
