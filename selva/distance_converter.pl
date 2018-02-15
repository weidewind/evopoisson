#!/usr/bin/perl

use Getopt::ArgvFile;
use Getopt::Long;

my $input;
my $output;
my $length;

GetOptions (
		'input=s' => \$input,
		'output=s' => \$output,
		'length=i' => \$length,
	);
	
	open FILE, "<$input" or die "cannot open $input\n";
	my $str = <FILE>;
	close FILE;
	$str =~ s!:(\d+)!":".$1/$length!ge;
	open OUT, ">$output" or die;
	print OUT $str;
	close OUT;
