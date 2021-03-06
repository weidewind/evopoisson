#!/usr/bin/perl

use Mutmap;

medtest();
#bintest();
funny();

sub medtest {
	my @array = (5,5,5,5);
	my %hash;
	$hash{0} = 0; # was skipped until 21.09.2016
	$hash{1} = 2; #until 21.09.2016 this hash had been turning into array (2,4,4,10), with 2 mutations at distance = 0
	$hash{2} = 4;
	$hash{3} = 4;
	$hash{4} = 10;
	print Mutmap::hist_mean(\@array,0.5)."\n";
	print Mutmap::hist_mean_for_hash(\%hash,0.5)."\n";
	print Mutmap::hist_mean(\@array)."\n";
	print Mutmap::hist_mean_for_hash(\%hash)."\n";
	print Mutmap::hist_median(\@array,0.5)."\n";
	print Mutmap::hist_median_for_hash(\%hash, 0.5)."\n";
	print Mutmap::hist_median(\@array)."\n";
	print Mutmap::hist_median_for_hash(\%hash)."\n";
	print Mutmap::hist_mean_for_hash(\%hash,2)."\n";
	print Mutmap::hist_median_for_hash(\%hash, 2)."\n";
}

sub funny {
	my %hash;
	#$hash{0} = 0; # was skipped until 21.09.2016
	$hash{1} = 5; #until 21.09.2016 this hash had been turning into array (2,4,4,10), with 2 mutations at distance = 0
	$hash{2} = 5;
	#$hash{2} = 4;
	#$hash{3} = 4;
	$hash{20} = 4;
	$hash{21} = 6;
	print Mutmap::hist_median_for_hash(\%hash,1)."\n";
	print Mutmap::hist_mean_for_hash(\%hash,1)."\n";
}
sub bintest {
	print Mutmap::bin(0, 1)."\t";
	print Mutmap::bin (1,1)."\t";
	print Mutmap::bin (1.2,1)."\t";
	print Mutmap::bin (1,2)."\t";
	print Mutmap::bin (1,0.5)."\t";
	print Mutmap::bin (2.5, 0.5)."\t";
	print Mutmap::bin (2, 10)."\t";
	print Mutmap::bin (0, 10)."\t";
}