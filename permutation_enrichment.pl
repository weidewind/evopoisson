#!/usr/bin/perl
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Mutmapnolim;
use Groups;
use Cwd qw(abs_path cwd getcwd);
use Sort::Rank  qw(rank_sort);
use Data::Dumper;
use List::Util  qw(sum min max);
use Enrichment qw(get_groups only_groups parse_real_stats diffvalues get_fake_dirs parse_fake_stats meanvalues);

$| = 1;

my $input;
my $protein;
my $data;
my $number;
#my $value = "pvalue_epistasis(mean)";

GetOptions (
		'data=s' => \$data,	
		'input=s' =>\$input,
		'protein=s' =>\$protein,
		'value=s' =>\$value,
		'number=i' =>\$number,
	);

my ($groups, $group_names) = get_groups($data,$protein);
my $all = $groups->[-1];
my ($g, $gn)  = only_groups($groups, $group_names);
my @sites = parse_real_stats($input,$protein);

print "$protein $value $input $data  \n";
my @values = ("pvalue_epistasis(mean)", "pvalue_epistasis(median)","pvalue_environment(mean)", "pvalue_environment(median)");
my %reals;
foreach my $v(@values){
	%{$reals{$v}} = meanvalues(\@sites,$g, $gn, $v);
}

print Dumper(\%reals);
my %enrichment_pvalues;
foreach my $i(1..$number){
	my @fakegroups;
	foreach my $group (@{$g}){
			my @fakegroup = permute($group, $all);
			push @fakegroups, \@fakegroup;
	}
#	print Dumper (\@fakegroups);
	foreach my $v(@values){
		%{$fakes{$v}} = meanvalues(\@sites,\@fakegroups, $gn, $v);
	}
#	print Dumper (\%fakes);
	foreach my $group (@{$gn}){
		foreach my $v(@values){
			if ($fakes{$v}{$group} <= $reals{$v}{$group}) {$enrichment_pvalues{$v}{$group} = $enrichment_pvalues{$v}{$group}+1;}
		}
	}
}
foreach my $group(@{$gn}){
	foreach my $v(@values){
		$enrichment_pvalues{$v}{$group} = $enrichment_pvalues{$v}{$group}/($number);
	}
}
 foreach my $gr(sort @{$gn}){
 	print $gr.",".$enrichment_pvalues{"pvalue_epistasis(mean)"}{$gr}.",".$enrichment_pvalues{"pvalue_environment(mean)"}{$gr}."\n";
 	print $gr.",".$enrichment_pvalues{"pvalue_epistasis(median)"}{$gr}.",".$enrichment_pvalues{"pvalue_environment(median)"}{$gr}."\n";
 }


sub permute{
	my $group = shift;
	my @all = @{$_[0]};
  	my @set;
  	for ( 1 .. scalar @{$group} ){
    	push @set, splice @all, rand @all, 1;
  	}
	return @set;
}
