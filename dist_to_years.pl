#!/usr/bin/perl 
use strict;
use warnings;
use Cwd qw(abs_path cwd getcwd);
use File::Spec;
use Bio::Phylo::IO;


my $tree = parse_tree("/export/home/popova/workspace/evopoisson/data/little_ksu/h3.l.r.newick");
my @terminals = @{ $tree->get_terminals };
open FILE, ">/export/home/popova/workspace/evopoisson/output/little_ksu/h3_dist_to_years";
foreach my $t (@terminals){
	my $name = $t->get_name;
	my @year = split(/_/, $name);
	my $dist = $t->calc_path_to_root();
	print FILE $year[1].",".$dist."\n";
}
close FILE;
my @internals = @{ $tree->get_internals };
open FILE, ">/export/home/popova/workspace/evopoisson/output/little_ksu/h3_node_dists";
foreach my $t (@internals){
	my $name = $t->get_name;
	my $dist = $t->calc_path_to_root();
	print FILE $name.",".$dist."\n";
}
close FILE;
sub parse_tree {
					my $tree_file = $_[0];
					open TREE, "<$tree_file" or die "Cannot open file ".$tree_file."\n";
					# get a newick string from some source
 					my $tree_string;
 					my $t;
 					while($t = <TREE>){
 						$t =~ s/\n$//;
 						$tree_string = $tree_string.$t;
 					}
 					 close TREE;

 					# Call class method parse from Bio::Phylo::IO
 					# note: newick parser returns 'Bio::Phylo::Forest'
                    # Call ->first to retrieve the first tree of the forest.
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $tree_string,
   					  -format => 'newick'
 					)->first;

 					return $tree;
	
} 