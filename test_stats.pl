#!/usr/bin/perl

#use Mutmapnolim qw(hist_median_for_hash hist_mean_for_hash hist_median hist_mean hist_to_array);
use List::Util qw(sum);
use strict;
use Data::Dumper;
use Mutmapnolim qw(tttplot);
use AgeingStat;
use BPStat;
use MeanStat;
use MedianStat;

#my $file = "/export/home/popova/workspace/evopoisson/testfiles/hashes";
#test($file); #  '1', '0.97752808988764'

my $file = "/export/home/popova/workspace/evopoisson/testfiles/allbins";
test($file); #  

sub test {
	my $file = shift;
	my $obshash;
	my $exphash;
	open(my $fh, '<', $file) or die "Can't read file '$file' [$!]\n";
	while (my $line = <$fh>) {
	    chomp $line;
	    my @fields = split(/,/, $line);
	    $obshash->{$fields[0]} = $fields[1];
	    $exphash->{$fields[0]} = $fields[2];
	}
#	my @plot = BPStat::tttplot($obshash, $exphash, 0.5);
#	foreach my $dot(@plot){
#		print $dot->[0].",".$dot->[1]."\n";
#	} 
	my $stat = AgeingStat->new("bp");
	$stat->computeStats({obshash=>$obshash, exphash=>$exphash});
	print "W is $stat\n";
}
