#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;
use MutMap::MutMapper;


my $mutmap = MutMap::MutMapper -> new("h1", "nsyn", "locally");
print_tree_for_sites($mutmap, 238);

# Procedure for printing trees
sub print_tree_for_sites {
	my $mutmap = $_[0];
	my @sites = @{$_[1]};
	for my $site(@sites){
		$mutmap -> print_static_tree_with_mutations($site);
	}
}