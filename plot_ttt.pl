#!/usr/bin/perl

#use Mutmapnolim qw(hist_median_for_hash hist_mean_for_hash hist_median hist_mean hist_to_array);
use List::Util qw(sum);
use strict;
use Data::Dumper;
use Getopt::Long;
use Getopt::ArgvFile;
use Mutmapnolim qw(tttplot);
use AgeingStat;
use BPStat;

my $file = "/export/home/popova/workspace/evopoisson/testfiles/hashes";
my $output;
my $norm;

GetOptions (	
		'input=s' => \$file,
		'output=s' => \$output,
		'norm' => \$norm,
	);
	
plot($file, $output, $norm); #  '1', '0.97752808988764'

sub plot {
	my $file = shift;
	my $output = shift;
	my $norm = shift;
	my ($obshash, $exphash) = readFile($file);
	my $stat = AgeingStat->new("bp");
	my $value = $stat->computeStats({obshash=>$obshash, exphash=>$exphash});
	my @plot = BPStat::tttplot($obshash, $exphash);
	if ($norm){
		@plot = BPStat::norm_plot(\@plot);
	}
	open FILE, ">$output" or die "Cannot open $output: $!\n";
	print "smth\n";
	foreach my $dot(@plot){
		print "smth else\n";
		print FILE $dot->[0].",".$dot->[1]."\n";
	} 
	close FILE;
	print "W is $value\n";
}

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