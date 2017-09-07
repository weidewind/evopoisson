#!/usr/bin/perl

use List::Util qw(sum);
use strict;
use Data::Dumper;
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Mutmapnolim qw(tttplot, concat, iterationFiles);

my $input = "/export/home/popova/workspace/evopoisson/output/h1test_no_mutnum_branch_skipnoskip_mutscontrolled_fixed/nsyn/maxpath_not_subtracted/h1/h1_for_enrichment_7";

GetOptions (	
		'input=s' => \$input,
	);

my %weeds = findWeeds($input);
printWeeds(\%weeds);

sub findWeeds {
	my $input = shift;
	if (-d $input){
		findWeedsInDir($input);
	}	
	else {findWeedsInFile($input);}
}

sub findWeedsInDir {
	my $dir = shift;
	my @files = Mutmapnolim::iterationFiles($dir);
	my $iterations;
	my %weeds;
	foreach my $file (@files){
			my %file_weeds = findWeedsInFile(File::Spec->catfile($dir,$file));
			$iterations += $file_weeds{iterations};
			foreach my $subtree (keys %{$file_weeds{weeds}}){
				$weeds{$subtree} += $file_weeds{weeds}{$subtree};
			}
	}
	return (iterations => $iterations, weeds => \%weeds);
}

sub findWeedsInFile {
	my $file = shift;
	open FILE, "<$file" or die "Cannot open $file: $!\n";
	my $iterations;
	my %weeds;
	while (<FILE>){
		if ($_ =~ /^>/){$iterations++;}
		if ($_ =~ /^s/){
			my %subtree = (split (/\s+/, $_));
			$weeds{Mutmapnolim::concat($subtree{site},$subtree{node})} += 1 if $subtree{muts} == 0;
		}
	}
	close FILE;
	return (iterations => $iterations, weeds => \%weeds);
}

sub printWeeds {
	my ($weeds) = @_;
	print $weeds->{iterations}."\n";
	foreach my $subtree(keys %{$weeds->{weeds}}){
		print $subtree." ".$weeds->{weeds}{$subtree}."\n";
	}
}
	
1;