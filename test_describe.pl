#!/usr/bin/perl
use List::Util qw(sum min);
use Cwd qw(abs_path cwd getcwd);
use File::Spec;

	my $outdir = "/export/home/popova/workspace/evopoisson/output/h1test_no_mutnum_branch_skipnoskip_mutscontrolled_fixed/nsyn/maxpath_not_subtracted/before";
 	opendir(DH, $outdir);
	my $prot = "h1";
	my @files = grep { /^${prot}_gulpselector_vector_boot_median_test_[0-9\.]+_single_sites/ } readdir(DH);
	closedir(DH);
	my @restr;
	foreach my $file(@files){
		$file =~ /.*_([0-9\.]+)_single_sites/;
		push @restr, $1;
	}
	@sorted = sort { $a <=> $b } @restr;
 	my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$sorted[0]."_single_sites");
 	print $file;
 