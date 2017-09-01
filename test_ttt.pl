#!/usr/bin/perl

#use Mutmapnolim qw(hist_median_for_hash hist_mean_for_hash hist_median hist_mean hist_to_array);
use List::Util qw(sum);
use strict;
use Data::Dumper;
use Mutmapnolim qw(tttplot);

my $file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/debugging/hashes";
my $obshash;
my $exphash;
open(my $fh, '<', $file) or die "Can't read file '$file' [$!]\n";
while (my $line = <$fh>) {
    chomp $line;
    my @fields = split(/,/, $line);
    $obshash->{$fields[0]} = $fields[1];
    $exphash->{$fields[0]} = $fields[2];
}

my $plot = tttplot($obshash, $exphash, 0.5);
print Dumper $plot;

sub test {
	
}