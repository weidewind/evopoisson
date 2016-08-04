#!/usr/bin/perl 
package Utils;
use strict;
use warnings;
use Bio::Phylo::IO;
use Bio::SeqIO;

our $VERSION = '0.01';


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

# read .fasta into a hash
sub parse_fasta {
	my $nodeseqs_file = shift;
	my %nodeseqs;
	my $seqio = Bio::SeqIO->new(-file => $nodeseqs_file, -format => "fasta");
	my $length;
	while ( my $seqobj = $seqio->next_seq ) {
		my $trimmed_id = (split(/\//, $seqobj->display_id))[0];
    	$nodeseqs{ $trimmed_id } = $seqobj->seq;
    	if (!$length){
    		$length = $seqobj->length();
    	}
	}
	return (\%nodeseqs, $length);
}

1;

=pod

=head1 SUPPORT

No support is available

=head1 AUTHOR

Copyright 2012 Anonymous.

=cut
