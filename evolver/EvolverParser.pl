#!/usr/bin/perl

use strict;
use Bio::Phylo::IO;
use EvolverParser;
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $dir;
my $tree;
GetOptions (	
		'dir=s' =>\$dir,
		'tree=s' =>\$tree,
	);

my $dirpath = File::Spec->catdir(getcwd(), $dir); 


my $ancs_file = = File::Spec->catfile($dirpath, "ancestral.txt");
my $tree_file = $tree; #"/export/home/popova/workspace/evopoisson/data/little_ksu/h1.l.r.newick";
my $strains_file = File::Spec->catfile($dirpath,"mc.paml");
my $outputdir = $dirpath;
my %hash = EvolverParser::evolverint_to_nodename($ancs_file, $tree_file);
foreach my $key (keys %hash){
	print $key." ".$hash{$key}."\n";
}
my @strains = EvolverParser::read_paml($strains_file);
my @ancs = EvolverParser::read_paml($ancs_file);

die "array sizes not equal!\n" if (scalar @strains != scalar @ancs);
for (my $i = 0; $i < scalar @strains; $i++){
	EvolverParser::print_fasta($strains[$i], $ancs[$i], \%hash, $outputdir."evolver_".$i.".fasta");
}