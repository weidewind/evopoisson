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

