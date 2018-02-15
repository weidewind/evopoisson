#!/usr/bin/perl

package EvolverParser;

use strict;
use Bio::Phylo::IO;



sub print_fasta {
	my %strains = %{$_[0]};
	my %ancs = %{$_[1]};
	my %evolverint_to_nodename = %{$_[2]};
	my $output = $_[3];
	open OUT, ">$output" or die "Cannot open $output for writing\n";
	foreach my $name (keys %strains){
		print OUT ">".$name."\n";
		print OUT $strains{$name}."\n";
	}
	foreach my $name (keys %ancs){
	print $name."\n";
	print substr($name,4)."\n";
		print OUT ">".$evolverint_to_nodename{substr($name,4)}."\n";
		print OUT $ancs{$name}."\n";
	}
}

sub read_paml {
	my $file = shift;
	open ST, "<$file" or die "Cannot open file ".$file."\n";
	my @sim_strains;
	my %fasta;
	for (1..2) {$_ = <ST>;}
	while(<ST>){
		next unless ($_ =~ /^[a-zA-Z]/);
		while ($_ && $_ =~ /^[a-zA-Z]/){
			my ($name, $seq) = read_pamlline($_);
			$fasta{$name} = $seq;
			$_ = <ST>;
		}
		push @sim_strains, {%fasta};
		my %fasta;
	}
	return @sim_strains;
		
}
sub read_pamlline{
	my $str = shift;
	my @splitter = split(/\s+/,$str);
	my $name = shift @splitter;
	my $seq;
	foreach my $block (@splitter){
		$seq = $seq.$block;
	}
	return ($name, $seq);
}

sub evolverint_to_nodename{
	my $ancs_file = shift;
	my $tree_file = shift;
	my %int_to_node;
	my %child_to_parent;
	 
	open TREE, "<$tree_file" or die "Cannot open file ".$tree_file."\n";
	my $tree_string = <TREE>;
	close TREE;
	my $tree = parse_tree($tree_file);
	
	open ANC, "<$ancs_file" or die "Cannot open file ".$ancs_file."\n";
	my $counter;
	my $outtree_string;
	my $branches;
	while(<ANC>){
		$counter++;
		$outtree_string = $_ if $counter == 3;
		$branches = $_ if $counter == 4;
	}
	close ANC;
	
	my @brs = split(/\s+/,$branches);
	foreach my $br(@brs){
		my @nodes = split(/\.{2}/, $br);
		$child_to_parent{$nodes[1]} = $nodes[0];
	}
	
	my @ints = ( $outtree_string =~ /[0-9]+/g );
	my @nodes = ( $tree_string =~ /STRAIN[0-9]+_[0-9]+/g );
	for (my $i = 0; $i < scalar @ints; $i++){
		$int_to_node{$ints[$i]} = $nodes[$i];
	}
	
	
	my %hash_of_nodes;
	my @nodes = $tree->get_nodes();
	foreach my $node(@nodes){
			my $name = $node ->get_name();
			$hash_of_nodes{$name} = \$node;
	}
	
	my %parents;
	foreach my $int (keys %int_to_node){
		my $nname = $int_to_node{$int};
		my $pint = $child_to_parent{$int};
		if ($pint){
			$parents{$pint} = 1;
			my $pname = ${$hash_of_nodes{$nname}}->get_parent()->get_name();
			$int_to_node{$pint} = $pname;
		}
	}
	
	while (%parents){
		my @parents = keys %parents;
		%parents = ();
		foreach my $int (@parents){
			my $nname = $int_to_node{$int};
			my $pint = $child_to_parent{$int};
			if ($pint){
				$parents{$pint} = 1;
				my $pname = ${$hash_of_nodes{$nname}}->get_parent()->get_name();
				$int_to_node{$pint} = $pname;
			}
		}
	}
	
	return %int_to_node;
}


sub parse_tree {
					my $tree_file = $_[0];
					open TREE, "<$tree_file" or die "Cannot open file ".$tree_file."\n";
 					my $tree_string;
 					my $t;
 					while($t = <TREE>){
 						$t =~ s/\n$//;
 						$tree_string = $tree_string.$t;
 					}
 					 close TREE;
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $tree_string,
   					  -format => 'newick'
 					)->first;
 					return $tree;
	
}
 
1;