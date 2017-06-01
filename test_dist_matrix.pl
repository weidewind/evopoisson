#!/usr/bin/perl

use List::Util qw(sum);
use Cwd qw(abs_path cwd getcwd);
use Bio::Phylo::IO;
use File::Spec;
use strict;

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

sub old_distance_matrix {
	    my $input_base = shift;
		my $file = File::Spec->catfile($input_base, "h3_distance_matrix.csv");
 		my %distance_hash;
	 	open CSV, "<$file" or die "Cannot open file $file\n";
	 	my $header = <CSV>;
	 	$header =~ s/[\s\n\r\t]+$//s;
	 	my @nodelables = split(',', $header);
	 	while(<CSV>){
				$_ =~ s/[\s\n\r\t]+$//s;
	 			my @dists = split(',', $_);
	 			my $node = $dists[0];
	 			for (my $i = 1; $i < scalar @dists; $i++){
	 				$distance_hash{$node}{$nodelables[$i]} = $dists[$i];
	 			}
	 	}
	 	close CSV;
	 	return %distance_hash;
}

sub new_distance_matrix {
	my $tree = shift;
	$tree->visit_breadth_first(
		-in   => sub{
			my $node=shift;
			if($node->is_root){
				$node->set_generic('time' => 0);
			}else{
				my $pnode=$node->get_parent;
				my $time=$pnode->get_generic('time');
				$time+=$node->get_branch_length;
				$node->set_generic('time' => $time);
			}
		}
	);	
}


sub old_distance {
	my $ancnode = shift;
	my $node = shift;
	my $distance_hash = shift;
	return $distance_hash->{$ancnode}->{$node};
}

sub new_distance {
	my $ancnode = shift;
	my $node = shift;
	my $hash_of_nodes = shift;
	$node = ${$hash_of_nodes->{$node}};
	$ancnode = ${$hash_of_nodes->{$ancnode}};
	return $node->get_generic('time') - $ancnode->get_generic('time');
}


my $input_base = File::Spec->catdir(getcwd(), "data", "ksu");
my $treefile = File::Spec->catfile($input_base, "h3.l.r.newick");
my $tree = parse_tree($treefile)  or die "No tree at $treefile";
my @nodes = $tree -> get_nodes;
my %hash_of_nodes;
foreach my $node(@nodes){
				my $name = $node ->get_name();
				$hash_of_nodes{$name} = \$node;
}


my %distance_hash = old_distance_matrix($input_base);
new_distance_matrix($tree);

sub test {
	my $ancnode = shift;
	my $node = shift;
	my $old = old_distance($ancnode, $node, \%distance_hash);
	print "old $old ";
	my $new = new_distance($ancnode, $node, \%hash_of_nodes);
	print " new $new \n";
	if ($old != $new){ print "Error!"}
}

#my $ancnode = "INTNODE3841";
#my $node = "INTNODE3840";
#test($ancnode, $node);

#my $ancnode = "INTNODE3841";
#my $node = "INTNODE3243";
#test($ancnode, $node);

my $ancnode = "Node_7476";
my $node = "Node_7466";
test($ancnode, $node);

