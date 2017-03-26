#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmap;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;


my $protein;
my $maxdepth = 50;
my $input = '';

GetOptions (	
		'protein=s' => \$protein,
		'input=s' => \$input,
		'maxdepth=i' => \$maxdepth,
	);
	
my $dirname = File::Spec->catdir(getcwd(), $input); 
my $filepath =  File::Spec->catfile($dirname,$protein."_sites");
open SITES, "<$filepath " or die "Cannot open $filepath : $!\n";
while (<SITES>){
	my @str = split(/,/, $_);
	push @{$hash{$str[1]}{$str[0]}}, \@str;
}
close SITES;
my @groups_and_names = Groups::only_groups($protein);
my @groups = @{$groups_and_names[0]};
my @names = @{$groups_and_names[1]};
my $output = File::Spec->catfile($dirname,$protein."_subtrees");
open GROUPSUBTREES, ">$output"; 
print scalar @{$groups};
for (my $i = 0; $i < scalar @groups; $i++){
	    print $names[$i]."\n";
		print GROUPSUBTREES $names[$i]."\n";
		for my $site (@{$groups[$i]}){
			if (exists $hash{$site}{$maxdepth}){
				my @strings = @{$hash{$site}{$maxdepth}};
				foreach my $str (@strings){
					print GROUPSUBTREES $site."\t".$str->[2]."\t".$str->[3]."\t".$str->[4]."\t".$str->[5]."\t".$str->[6]."\t".$str->[7]."\t".$str->[8]."\n";
				}
			}
		}
}
close GROUPSUBTREES;


