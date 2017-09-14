#!/usr/bin/perl
package Weeds;
use List::Util qw(sum);
use strict;
use Data::Dumper;
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Textbits qw(concat cleave iterationFiles);



sub new {
	my $class = shift;
	my $input = shift;
	my $weeds;
	if (-d $input){
		$weeds = findWeedsInDir($input);
	}	
	else {$weeds = findWeedsInFile($input);}
	$weeds->{fails_threshold} = 0;
	return bless $weeds, $class;	
}

sub findWeedsInDir {
	my $dir = shift;
	my @files = Textbits::iterationFiles($dir);
	my $iterations;
	my %weeds;
	foreach my $file (@files){
			my $file_weeds = findWeedsInFile(File::Spec->catfile($dir,$file));
			$iterations += $file_weeds->{iterations};
			foreach my $subtree (keys %{$file_weeds->{weeds}}){
				$weeds{$subtree} += $file_weeds->{weeds}{$subtree};
			}
	}
	return {iterations => $iterations, weeds => \%weeds};
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
			my $sn = Textbits::concat($subtree{site},$subtree{node});
			$weeds{$sn} += 1 if $subtree{muts} == 0;
		}
	}
	close FILE;
	return {iterations => $iterations, weeds => \%weeds};
}

sub printWeeds{
	my $weeds = shift;
	my $outfile = shift;
   	open(OUTPUT, '>', $outfile) or die "Cannot open $outfile for writing: $!\n";
	print OUTPUT "Total number of iterations: #".$weeds->{iterations}."\n";
	print OUTPUT "Fails threshold: #".$weeds->{fails_threshold}."\n";
	print OUTPUT "site:node,fails_number\n";
	foreach my $subtree(keys %{$weeds->{weeds}}){
		print OUTPUT $subtree.",".$weeds->{weeds}{$subtree}."\n";
	}
	close OUTPUT;
}

sub readWeeds{
	my $class = shift;
	my $weeds_file = shift;
	my $weeds;
	open WEEDS, "<$weeds_file" or die "Cannot open weeds file $weeds_file:$!\n";
	my $str = <WEEDS>;
	chomp ($weeds->{iterations} = (split(/#/, $str))[1]);
	$str = <WEEDS>;
	chomp ($weeds->{fails_threshold} =  (split(/#/, $str))[1]);
	$str = <WEEDS>;
	while(<WEEDS>){
		chomp;
		my ($subtree, $fails) = split(/,/, $_);
		$weeds->{weeds}{$subtree} =  $fails;
	}
	close WEEDS;
	return bless $weeds, $class;
}

sub worstWeeds{
	my $weeds = shift;
	my ($args) = @_;
	my $failsthr = $args->{fails_threshold};
	my $its =  $weeds->{iterations};
	my $worst;
	$worst->{fails_threshold} = $failsthr;
	$worst->{iterations} = $its;
	foreach my $subtree(keys %{$weeds->{weeds}}){
		my $fails = $weeds->{weeds}{$subtree};
		if (($its-$fails)/$its < $failsthr){
			$worst->{weeds}{$subtree} = $fails;
		}
	}
	return bless $worst, ref $weeds;
}
	
1;