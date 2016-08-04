#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

use strict;

my %hash = (a => 1, b => 2, c=> 3);
change1(\%hash);
print $hash{a};
change2(\%hash);
print $hash{a};
change3(\%hash);
print $hash{a};

sub change1 {
	my %hash = %{$_[0]};
	$hash{a} = 0;
}

sub change2 {
	$_[0] -> {a} = 7;
}

sub change3 {
	my $hashref = $_[0];
	$hashref -> {a} = 9;
}


