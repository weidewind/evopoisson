#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use MutMap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;



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

my $mutmap = MutMap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state, fromfile => true});
my @groups_and_names = $mutmap -> get_groups_and_names_for_protein(); 

## for concat_and_divide_simult you need a mutmap produced from realdata, therefore fromfile => true

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