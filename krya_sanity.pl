#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Getopt::ArgvFile;


our @h1_leading_kr = (4,11,13,16,52,60,73,74,86,88,90,91,97,99,111,113,128,144,151,156,157,162,169,170,171,172,173,178,182,184,199,202,203,205,206,207,209,210,232,240,261,268,269,283,287,289,290,293,326,361,415,488,489);
our @h1_trailing_kr = (3,6,7,11,52,53,64,89,91,99,101,111,129,137,142,144,148,150,156,157,158,162,165,169,172,178,179,185,186,195,197,199,200,201,203,207,231,236,243,251,253,269,274,278,289,290,293,299,324,331,361,389,390,398,415,422,455,467,470,489,510,514,515,526,562);
my @h3_leading_kr = (11,19,41,66,73,78,91,99,122,137,140,147,149,151,153,156,159,160,161,171,172,174,175,176,179,188,202,205,206,208,209,213,218,233,235,238,242,243,277,278,291,292,294,363,377,391,402,466,468,545);
my @h3_trailing_kr = (3,10,11,13,16,43,49,61,63,64,65,66,67,69,72,73,91,101,104,108,110,112,122,128,144,147,153,154,156,158,166,173,175,176,179,180,187,188,189,190,208,209,210,215,217,218,223,225,230,232,236,238,239,241,242,245,260,264,276,285,289,292,295,307,339,342,377,391,401,402,434,442,466,484,505,516,538,546,563);
our @n2_leading_kr = (18,20,23,30,52,93,143,150,194,197,199,208,216,220,221,249,265,307,308,310,313,328,336,339,344,346,368,369,370,372,381,385,387,390,432);
our @n2_trailing_kr = (2,4,5,9,27,30,40,44,45,50,56,65,77,82,83,120,127,147,148,149,151,155,210,216,220,238,248,251,258,263,265,269,302,307,309,310,312,328,329,334,335,338,339,342,347,372,386,392,400,402,403,414,416,432,433,434,455,464);
our @n1_leading_kr = (15,17,23,34,45,64,70,78,105,173,200,214,222,234,248,249,250,254,270,274,275,287,329,332,336,339,344,352,354,367,369,382,390,396,418,427,430,434,451);
our @n1_trailing_kr = (15,17,21,23,38,39,40,42,45,47,48,52,57,67,68,70,73,77,81,82,83,93,100,101,114,130,147,149,188,200,249,254,259,262,264,267,270,273,275,329,331,340,346,352,364,366,367,390,416,418,419,427,435,452,453,455,462);

my %hash;
$hash{"h1"}{"leading"} = \@h1_leading_kr;
$hash{"h1"}{"trailing"} = \@h1_trailing_kr;
$hash{"h3"}{"leading"} = \@h3_leading_kr;
$hash{"h3"}{"trailing"} = \@h3_trailing_kr;
$hash{"n1"}{"leading"} = \@n1_leading_kr;
$hash{"n1"}{"trailing"} = \@n1_trailing_kr;
$hash{"n2"}{"leading"} = \@n2_leading_kr;
$hash{"n2"}{"trailing"} = \@n2_trailing_kr;

foreach my $protein (keys %hash){
	foreach my $gr ("trailing", "leading"){
		foreach my $state ("nsyn", "syn"){
			my $args = {protein => $protein, state => $state}; 
			my $mutmap = Mutmap->new($args);
			foreach my $site(@{$hash{$protein}{$gr}}){
				print "$site $protein $state $gr\n";
				unless (exists $mutmap->{static_nodes_with_sub}{$site}) {
					print "no muts at $site $protein $state $gr !\n";
				}
			
			}
		}
	}
}