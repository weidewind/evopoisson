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
use Statistics::Descriptive;

$| = 1;

my $input;
my $protein;
my $data;
my $fakesdir;
my $tail;
#my $value = "pvalue_epistasis(mean)";

GetOptions (
		'data=s' => \$data,	
		'input=s' =>\$input,
		'protein=s' =>\$protein,
		'value=s' =>\$value,
		'fakesdir=s' => \$fakesdir,
		'tail=s' => \$tail,
	);

my ($groups, $group_names) = get_groups($data,$protein);
my ($g, $gn)  = only_groups($groups, $group_names);
my @sites = parse_real_stats($input,$protein);

print "$protein $value $input $data $fakesdir \n";
my @values = ("pvalue_epistasis(mean)", "pvalue_epistasis(median)","pvalue_environment(mean)", "pvalue_environment(median)");
my %reals;
foreach my $v(@values){
	%{$reals{$v}} = diffvalues(\@sites,$groups, $group_names, $v);
}

my @fakedirs = get_fake_dirs($fakesdir, $tail);
my %enrichment_pvalues;
my %fakeholder;
foreach my $dir(@fakedirs){
	my @fakesites = parse_fake_stats($dir,$protein);
	my %fakes;
	foreach my $v(@values){
		%{$fakes{$v}} = diffvalues(\@fakesites,$groups, $group_names, $v);
	}

	foreach my $group (@{$gn}){
		foreach my $v(@values){
			push @{$fakeholder{$v}{$group}}, $fakes{$v}{$group};
			if ($fakes{$v}{$group} <= $reals{$v}{$group}) {$enrichment_pvalues{$v}{$group} = $enrichment_pvalues{$v}{$group}+1;}
		}
	}
}
foreach my $group(@{$gn}){
	foreach my $v(@values){
		$enrichment_pvalues{$v}{$group} = $enrichment_pvalues{$v}{$group}/(scalar @fakedirs);
	}
}
 foreach my $gr(sort @{$gn}){
 	print $gr.",".$enrichment_pvalues{"pvalue_epistasis(mean)"}{$gr}.",".$enrichment_pvalues{"pvalue_environment(mean)"}{$gr}."\n";
 	print $gr.",".$enrichment_pvalues{"pvalue_epistasis(median)"}{$gr}.",".$enrichment_pvalues{"pvalue_environment(median)"}{$gr}."\n";
 }
	
print "group,value,real,fake_mean,fake_std,lower95,upper95\n";	
 foreach my $gr(sort @{$gn}){ 
 	foreach my $v(@values){
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(\@{$fakeholder{$v}{$gr}});
		my $upper = $stat->percentile(95);
		my $lower = $stat->percentile(5);
		print $gr.",".$v.",".$reals{$v}{$gr}.",".$stat->mean().",".$stat->standard_deviation().",".$lower.",".$upper."\n";
 	}
 }



