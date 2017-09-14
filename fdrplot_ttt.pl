#!/usr/bin/perl

#use Mutmapnolim qw(hist_median_for_hash hist_mean_for_hash hist_median hist_mean hist_to_array);
use List::Util qw(sum);
use strict;
use Data::Dumper;
use Getopt::Long;
use Getopt::ArgvFile;
use FileTTT;
use Mutmapnolim qw(tttplot);
use AgeingStat;
use BPStat;

my $dir = "/export/home/popova/workspace/evopoisson/output/fakes_h1_no_mutnum_control_branch_controlled_fixed";
my $innerpath = "nsyn/maxpath_not_subtracted/skip_stoppers";
my $file ="h1_gulpselector_vector_boot_median_test_0_all";
my $output = "output/tttplots/testplot";
my $norm;

GetOptions (	
		'dir=s' => \$dir,
		'innerpath=s' => \$innerpath,
		'file=s' => \$file,
		'output=s' => \$output,
		'norm' => \$norm,
	);
	
folderplot($dir, $innerpath, $file, $output, $norm); #  '1', '0.97752808988764'

sub folderplot {
	my $dir = shift;
	my $innerpath = shift;
	my $filename = shift;
	my $output = shift;
	my $norm = shift;
	
	opendir(DH, $dir);
	my @fakedirs = readdir(DH);
	closedir(DH);
	die "No files or dirs found in $input!\n" unless scalar @fakedirs > 0;
	my %files;
	foreach my $di(sort @fakedirs){
		my $dipath =  File::Spec->catdir($dirname,$di);
		if (-d $dipath && $di =~ /([0-9]+)_fake/){
			my $fakeno = $1;
			$files{$fakeno} = File::Spec->catfile($dipath,$innerpath,$filename);
		}
	}
	open OUT, ">$output" or die "Cannot open $output: $!\n";
	foreach my $fakeno(keys %files){
		my ($obshash, $exphash) = readFile($files{$fakeno});
		my @plot = BPStat::tttplot($obshash, $exphash);
		if ($norm){
			@plot = BPStat::norm_plot(\@plot);
		}
		foreach my $dot(@plot){
			print OUT $fakeno.",".$dot->[0].",".$dot->[1]."\n";
		} 
	}
	close OUT;
}