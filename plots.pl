#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use MutMap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);



my $protein;
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $verbose;
GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'verbose'  => \$verbose,
	);

my $mutmap = MutMap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state});

## prints files for drawing the chosen version of plots: mutations_density_in_circles(distance from ancestral mutation), points correspond to mutations
## There could be more than one mutation at the given distance, so number of points on the plot can be less than the number of mutations in the subtree under analysis 
$mutmap->egor_smart_site_entrenchment($verbose);



#prepare_real_data("h1",0,0,0,"locally");
#prepare_real_data("h3",0,0,0,"locally");
#prepare_real_data("n1",0,0,0,"locally");
#prepare_real_data("n2",0,0,0,"locally");
#my $realdata = lock_retrieve ("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1_realdata") or die "Cannot retrieve real_data";
#my $norm100 = compute_norm(100, \@h1_antigenic);
#print "\nnorm 100 ".$norm100;
#my $norm150 = compute_norm(150, \@h1_antigenic);
#print "\nnorm 150 ".$norm150;
#my $obshash50 = select_obshash(150, \@h1_antigenic);

## Procedure for testing gulp_iterations locally
#my $prot = $ARGV[0];
#my $tag = "locallength_test";
#my $iterations = 3; 
#$static_protein = $prot;
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata") or die "Cannot retrieve real_data";
#print "Starting gulp $tag of $iterations iterations for protein $prot..\n";
#iterations_gulp ($prot, $iterations, 0, 1, $tag, 1);
#print "Finished gulp $tag of $iterations iterations for protein $prot\n";
#test_gulp("/cygdrive/c/Users/weidewind/Documents/CMD/misc_scripts/h3_for_enrichment_locallength_test", 150);
###




## Procedure for launching iterations on server
#my $prot = $ARGV[0];
#my $tag = $ARGV[1];
#my $iterations = 500; 
#$static_protein = $prot;
##my $realdata = lock_retrieve ("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_realdata_maxpath_not_subtracted") or die "Cannot retrieve real_data";
#my $realdata = lock_retrieve ("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_realdata_maxpath_not_subtracted") or die "Cannot retrieve real_data";
#print "Starting gulp $tag of $iterations iterations for protein $prot..\n";
#iterations_gulp ($prot, $iterations, 0, 0, $tag);
#print "Finished gulp $tag of $iterations iterations for protein $prot\n";
###

# Test procedure for checking concat_and_divide
#$static_protein = "h1";  
#
#my $tag = "tag1";
#print "stat subtree info check ";
#print scalar keys %static_subtree_info;
#print "\n";
#	my $realdata = lock_retrieve ("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_realdata") or die "Cannot retrieve real_data";
#
#iterations_gulp ($prot, $iterations, 0, 1, $tag); # in fact, it does not read restriction parameter (third)
#$tag = "tag2";
#iterations_gulp ($prot, $iterations, 0, 1, $tag);
#$tag = "tag3";
#iterations_gulp ($prot, $iterations, 0, 1, $tag);
#my @norms = (10,10,10);
#concat_and_divide ($prot, 3, "tag", \@norms);
#count_pvalues($prot,"tag");
##

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
#my @groups;
#my @names;
#if ($prot eq "h1"){
#	@groups = (\@h1_increased_binding, \@h1_antigenic, \@h1_pocket_closest, \@h1_surface, \@h1_internal, \@h1_host_shift_001, \@h1_leading_kr, \@h1_trailing_kr, \@h1_antigenic_ren);
#	@names = ("increased_binding", "antigenic", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr", "antigenic_ren");
#}
#elsif ($prot eq "h3"){
#	@groups = (\@h3_shih_epitopes, \@h3_antigenic, \@h3_antigenic_koel, \@h3_pocket_closest, \@h3_surface, \@h3_internal, \@h3_host_shift_001, \@h3_leading_kr, \@h3_trailing_kr);
#	@names = ("shih_epitopes", "antigenic", "antigenic_koel", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr");
#}
#elsif ($prot eq "n1"){
#	@groups = (\@n1_epitopes, \@n1_wan_epitopes, \@n1_pocket_closest, \@n1_surface, \@n1_internal, \@n1_host_shift_001, \@n1_leading_kr, \@n1_trailing_kr);
#	@names = ("epitopes", "wan_epitopes", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr");
#}
#elsif ($prot eq "n2"){
#	@groups = (\@n2_epitopes, \@n2_pocket_closest, \@n2_surface, \@n2_internal, \@n2_host_shift_001, \@n2_leading_kr, \@n2_trailing_kr, \@n2_decreasing, \@n2_increasing);
#	@names = ("epitopes", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr", "decreasing", "increasing");
#}
#my @groups_and_names = prepare_groups_and_names(\@groups, \@names);
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



## 26.02 Number of nodes in analysis
#print "\nh1------\n";
#my $prot = "h1";
#my $restriction = 50;
#my @all = (0..565);
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @groups = (\@all, \@h1_increased_binding, \@h1_antigenic_ren, \@h1_epitopes, \@h1_antigenic, \@h1_pocket_closest, \@h1_surface, \@h1_internal, \@h1_host_shift_001, \@h1_leading_kr, \@h1_trailing_kr);
#my @names = ("all", "increased_binding", "antigenic_ren","epitopes", "antigenic","pocket_closest","surface","internal","host_shift_001","leading", "trailing");
#print_nodes_in_analysis ($prot, $restriction, \@groups, \@names, 1);
#$realdata = 0;
#print "\nh3------\n";
#my $prot = "h3";
#my $restriction = 50;
#my @all = (0..565);
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @groups = (\@all, \@h3_antigenic, \@h3_antigenic_koel, \@h3_antigenic_smith, \@h3_antigenic_steinbruck, \@h3_epitopes, \@h3_shih_epitopes, \@h3_pocket_closest, \@h3_surface, \@h3_internal, \@h3_host_shift_001, \@h3_leading_kr, \@h3_trailing_kr);
#my @names = ("all", "antigenic", "antigenic_koel", "antigenic_smith", "antigenic_steinbruck", "epitopes", "shih_epitopes","pocket_closest","surface","internal","host_shift_001","leading", "trailing");
#print_nodes_in_analysis ($prot, $restriction, \@groups, \@names, 1);
#$realdata = 0;
#print "\nn2------\n";
#my $prot = "n2";
#my $restriction = 50;
#my @all = (0..565);
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @groups = (\@all, \@n2_epitopes, \@n2_pocket_closest, \@n2_surface, \@n2_internal, \@n2_host_shift_001, \@n2_leading_kr, \@n2_trailing_kr, \@n2_decreasing, \@n2_increasing);
#my @names = ("all", "epitopes","pocket_closest","surface","internal","host_shift_001","leading", "trailing", "decreasing", "increasing");
#print_nodes_in_analysis ($prot, $restriction, \@groups, \@names, 1);
#$realdata = 0;
#print "\nn1------\n";
#my $prot = "n1";
#my $restriction = 50;
#my @all = (0..565);
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @groups = (\@all, \@n1_wan_epitopes, \@n1_epitopes, \@n1_pocket_closest, \@n1_surface, \@n1_internal, \@n1_host_shift_001, \@n1_leading_kr, \@n1_trailing_kr);
#my @names = ("all", "wan_epitopes","epitopes", "pocket_closest","surface","internal","host_shift_001","leading", "trailing");
#print_nodes_in_analysis ($prot, $restriction, \@groups, \@names, 1);



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


# 18.02 syn
#my $prot = "h3";
#my $tag = "synall";
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/syndata/".$prot."_syndata");
#my @maxdepths = (50, 100, 150);
#my @groups = (\@h3_epitopes,\@h3_antigenic_smith,\@h3_antigenic_steinbruck, \@h3_antigenic_koel, \@h3_pocket_closest, \@h3_surface, \@h3_internal, \@h3_host_shift_001, \@h3_leading_kr, \@h3_trailing_kr);
#my @names = ("epitopes", "antigenic_smith", "antigenic_steinbruck","antigenic_koel","pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr", "antigenic_ren");
##my @groups = (\@h1_epitopes);
##my @names = ("epitopes");
#my @groups_and_names = prepare_groups_and_names(\@groups, \@names);
#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}, "syn");
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});

# 18.02 syn
#my $prot = "n2";
#my $tag = "synall";
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/syndata/".$prot."_syndata");
#my @maxdepths = (50, 100, 150);
#my @groups = (\@n2_internal);
#my @names = ("internal");
#my @groups_and_names = prepare_groups_and_names(\@groups, \@names);
#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]}, "syn");
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});

# 16.02 leading to trailing enrichment
#my $prot = "n2";
#my $tag = "lt";
#my $realdata = lock_retrieve ("/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_realdata");
#my @maxdepths = (50, 100, 150);
#my @groups_and_names = prepare_lt_groups_and_names($prot);
##print scalar @{$groups_and_names[0]};
##print "\n";
##print scalar @{$groups_and_names[1]};
##print "\n";
##foreach my $n(@{$groups_and_names[0]}){
##	foreach my $t(@{$n}){
##		print $t."\n";
##	}
##}
##print scalar @{$groups_and_names[1]};
##print "\n";
##foreach my $n(@{$groups_and_names[1]}){
##	print $n."\n";
##}
#concat_and_divide_simult ($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});
#count_pvalues($prot, $tag, \@maxdepths, \@{$groups_and_names[0]}, \@{$groups_and_names[1]});




# 18.02 syn control
#prepare_syn_data("n2", 0,0,1);

