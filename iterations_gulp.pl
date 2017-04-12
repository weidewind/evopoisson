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
my $tag;		# there will be several gulp files in one folder, therefore a tag for each gulp must be specified
my $iterations = 500;
my $subtract_tallest = '0';
my $verbose;
my $memusage;
my $no_neighbour_changing;
my $no_leaves;
my $restriction;
my $faketag;

my $lifetime_restr;
my $onestrip;
my $shuffler_type = "exp";
my $debugmode;
my $poisson;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'tag=s' => \$tag,
		'iterations=i' =>\$iterations,
		'subtract_tallest=i' => \$subtract_tallest,
		'verbose'  => \$verbose,
		'memusage'  => \$memusage, # should i print memusage in file?
		'no_neighbour_changing' => \$no_neighbour_changing,
		'no_leaves' => \$no_leaves,
		'restriction=i' => \$restriction,
		'faketag=s' => \$faketag,
		'lifetime_restr' => \$lifetime_restr,
		'onestrip' => \$onestrip,
		'shuffler_type=s' => \$shuffler_type,
		'debugmode' => \$debugmode,
		'distrpoisson' =>\$poisson,
	);


## Procedure for launching a gulp of iterations
unless ($subtract_tallest == 0 || $subtract_tallest == 1) {die "--subtract_tallest must be either 0 or 1\n";}
unless (defined $tag){die "There will be several gulp files in one folder, therefore a --tag for each gulp must be specified (and it should be an integer)";}
my $args = {bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, subtract_tallest => $subtract_tallest, no_neighbour_changing => $no_neighbour_changing, no_leaves => $no_leaves, fromfile => 1, faketag => $faketag}; 
unless  (Mutmapnolim::realdata_exists($args)) { die "No such realdata!"; }
print "realdata restriction is ".Mutmapnolim::check_realdata_restriction($args)."\n";

if ($verbose) { print "Starting gulp $tag of $iterations iterations for protein $protein..\n"; }
## for launching iterations you need a mutmap produced from realdata, therefore fromfile => true
my $mutmap = Mutmapnolim->new($args);
$mutmap-> iterations_gulp ($iterations, $tag, $verbose, $memusage, $restriction, $lifetime_restr, $onestrip, $shuffler_type, $debugmode, $poisson);
if ($verbose) { print "Finished gulp $tag of $iterations iterations for protein $protein\n"; }
###



## 25.01 Procedure for obtaining p-values


#my $prot = "h3";
#my $tag = "26_02";
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @maxdepths = (50, 100, 150);
#my @groups = (\@h3_shih_epitopes);
#my @names = ("shih_epitopes");
##my @groups = (\@h3_antigenic, \@h3_antigenic_koel, \@h3_pocket_closest, \@h3_surface, \@h3_internal, \@h3_host_shift_001, \@h3_leading_kr, \@h3_trailing_kr);
##my @names = ("antigenic", "antigenic_koel", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr");
#my @groups_and_names = prepare_groups_and_names(\@groups, \@names);
#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});


#my $prot = "n1";
#my $tag = "19_02";
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @maxdepths = (50, 100, 150);
##my @groups = (\@n1_wan_epitopes, \@n1_pocket_closest, \@n1_surface, \@n1_internal, \@n1_host_shift_001, \@n1_leading_kr, \@n1_trailing_kr);
##my @names = ("wan_epitopes", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr");
#my @groups = (\@n1_epitopes);
#my @names = ("epitopes");
#my @groups_and_names = prepare_groups_and_names(\@groups, \@names);
#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});

#my $prot = "n2";
#my $tag = "15_02";
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @maxdepths = (50, 100, 150);
#my @groups = (\@n2_epitopes, \@n2_pocket_closest, \@n2_surface, \@n2_internal, \@n2_host_shift_001, \@n2_leading_kr, \@n2_trailing_kr, \@n2_decreasing, \@n2_increasing);
#my @names = ("epitopes", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr", "decreasing", "increasing");
#my @groups_and_names = prepare_groups_and_names(\@groups, \@names);
#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});


# 18.02 syn
#my $prot = "h1";
#my $tag = "synall";
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/syndata/".$prot."_syndata");
#my @maxdepths = (50, 100, 150);
#my @groups = (\@h1_epitopes, \@h1_antigenic, \@h1_pocket_closest, \@h1_surface, \@h1_internal, \@h1_host_shift_001, \@h1_leading_kr, \@h1_trailing_kr, \@h1_antigenic_ren);
#my @names = ("epitopes", "antigenic", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr", "antigenic_ren");
##my @groups = (\@h1_epitopes);
##my @names = ("epitopes");
#my @groups_and_names = prepare_groups_and_names(\@groups, \@names);
#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}, "syn");
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});