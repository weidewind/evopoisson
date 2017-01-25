#!/usr/bin/perl

use strict;
use Bio::Phylo::IO;
use shuffle_muts_on_tree qw(shuffle_mutations_on_tree prepare_shuffler StripConstrains Shuffler);

my $tree = parse_tree("/export/home/popova/workspace/evopoisson/data/h1.l.r.newick");
my @nodes = $tree -> get_nodes;
my %hash_of_nodes;	
foreach my $node(@nodes){
	#if ($node->is_root()) {next;}
	my $name = $node ->get_name();
	$hash_of_nodes{$name} = $node;
}

my $ancestor_nodes = {
	INTNODE2429 => 1,
	INTNODE2406 => 1,
};

my $rh_constrains = {
	INTNODE2429 => {
		171 => StripConstrains->new(number_of_mutations => 6,lifetime => 230, stoppers => [$hash_of_nodes{"INTNODE2285"},$hash_of_nodes{"INTNODE1538"}]),  # stoppers are mocked
		184 => StripConstrains->new(number_of_mutations => 5,lifetime => 230, stoppers => [$hash_of_nodes{"INTNODE2285"},$hash_of_nodes{"INTNODE1354"}]),
	},
	INTNODE2406 => {
		202 => StripConstrains->new(number_of_mutations => 7,lifetime => 106, stoppers => [$hash_of_nodes{"INTNODE2394"}]),
		418 => StripConstrains->new(number_of_mutations => 6,lifetime => 143, stoppers => [$hash_of_nodes{"INTNODE2301"}]),
		205 => StripConstrains->new(number_of_mutations => 9,lifetime => 123, stoppers => [])
	}
	
};

#my $rh_constrains = $self->get_constraints($ancestor_nodes); #todo
	# $rh_constrains->{$name}->{$site}->StripConstrains
	#struct StripConstrains =>{
	#number_of_mutations => '$',
	#lifetime => '$',
	#stoppers => '@'
	#}
	my $shuffler = shuffle_muts_on_tree::prepare_shuffler($tree, $rh_constrains); #todo	
	my $rh_out_subtree = shuffle_muts_on_tree::shuffle_mutations_on_tree($shuffler); #todo
	
	
	foreach my $name (keys %{$rh_out_subtree}){
		print $name."\n";
		#foreach my $s (@{$shuffler->sites_hash($name)}){
		#	print "$s"."\t";
		#}
		print "\n";
		foreach my $site (keys %{$rh_out_subtree->{$name}}){
			print "\t".$site."\n";
			print "\t\t";
			foreach my $exitnode (@{$rh_out_subtree->{$name}->{$site}}){
				print $exitnode."\t";
			}
			print "\n";
		}
	}
	# $rh_out_subtree->{$name}->{$site}= массив уходов #todo
	
	sub parse_tree {
					my $tree_file = $_[0];
					open TREE, "< $tree_file" or die "Cannot open file ".$tree_file."\n";
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