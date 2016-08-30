#!/usr/bin/perl

use lib getcwd(); # adds working directory to @INC
use MutMap;
use Getopt::Long;
use Groups;



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $tag;
my $iterations = 500;
my $verbose;

GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'tag=s' => \$tag,
		'iterations=i' =>\$iterations,
		'verbose'  => \$verbose,
	);


## Procedure for launching a gulp of iterations

 
if ($verbose) { print "Starting gulp $tag of $iterations iterations for protein $prot..\n"; }
## for launching iterations you need a mutmap produced from realdata, therefore fromfile => true
my $mutmap = MutMap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, fromfile => true});
$mutmap-> iterations_gulp ($iterations, $tag);
if ($verbose) { print "Finished gulp $tag of $iterations iterations for protein $prot\n"; }
###



## 25.01 Procedure for obtaining p-values
#my $prot = "n1";
#my $syn = 0;
#my $subtract_maxpath = 1;
#my $tag = "median_vs_mean"."_".syn_tag($syn)."_".maxpath_tag($subtract_maxpath);
#
##my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
### For this line to work with local data you will have to change the structure
#my $dir = File::Spec->catdir(getcwd(),syn_tag($syn), maxpath_tag($subtract_maxpath));
#my $realdatapath = File::Spec->catfile($dir, $prot."_realdata");
#my $realdata = lock_retrieve ($realdatapath);
#my @maxdepths = (50, 100, 150);

my @groups_and_names = $mutmap -> get_groups_and_names_for_protein(); 

#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}, $syn, $subtract_maxpath);
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}, $dir);

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