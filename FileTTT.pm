#!/usr/bin/perl

#use Mutmapnolim qw(hist_median_for_hash hist_mean_for_hash hist_median hist_mean hist_to_array);
use List::Util qw(sum);
use strict;

sub readFile {
	my $file = shift;
	my $obshash;
	my $exphash;
	open IN, "<$file" or die "Can't read file '$file' [$!]\n";
	my $str = <IN>;
	close IN;
	if ($str =~ /^bin/){
		my ($obshash, $exphash) = readPvalFile($file);
	}
	elsif(grep (/,/, $str)){
		my ($obshash, $exphash) = readCSV($file);
	}	
	
}

sub readPvalFile {
	my $file = shift;
	my $obshash;
	my $exphash;
	open IN, "<$file" or die "Can't read file '$file' [$!]\n";
	while (<IN>) {
		if ($_ =~ /^$/){last;}
		if ($_ =~ /!^[0-9]/){next;}
	    chomp $_;
	    print $_."\n";
	    my @fields = split(/[\t\s]+/, $_);
	    $obshash->{$fields[0]} = $fields[1];
	    $exphash->{$fields[0]} = $fields[2];
	}
	close IN;
	return ($obshash, $exphash);
}

sub readCSV {
	my $file = shift;
	my $obshash;
	my $exphash;
	open IN, "<$file" or die "Can't read file '$file' [$!]\n";
	while (<IN>) {
	    chomp $_;
	    my @fields = split(/,/, $_);
	    $obshash->{$fields[0]} = $fields[1];
	    $exphash->{$fields[0]} = $fields[2];
	}
	close IN;
	return ($obshash, $exphash);
}

1;