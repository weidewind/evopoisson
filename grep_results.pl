#!/usr/bin/perl
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $input;
GetOptions (	
		'input=s' =>\$input,
	);

my $dirname = File::Spec->catdir(getcwd(), $input); 
opendir(DH, $dirname);
my @files = readdir(DH);
closedir(DH);
if (scalar @files == 0){
	die "No files found in $input!\n";
}

foreach my $filename(sort @files){
	my $filepath =  File::Spec->catfile($dirname,$filename);
	next if (-d $$filepath);
	next unless ($filename =~ /(.*)_gulpselector_vector_boot_median_test_([0-9]+)_(.*)/);
	my $prot = $1;
	my $depth = $2;
	my $ending = $3;
	open FILE, "<$filepath" or die "Cannot open $filename: $!\n";
	if ($ending =~ /(.*)_complement$/ || $ending =~ /(all)$/){
		my $group = $1;
		my $groupfile = File::Spec->catfile($dirname, $prot."_groups");
		open GROUPS, ">>$groupfile" or die "Cannot open $groupfile: $!\n";
		while(<FILE>){
			if ($_ =~ /^me/){
				print GROUPS $group."\t".$depth."\t".$_;
			}
		}
		close GROUPS;
		 
	}
	elsif ($ending =~ /.*single_sites$/){
		my $sitesfile = File::Spec->catfile($dirname, $prot."_sites");
		open SITES, ">>$sitesfile" or die "Cannot open $sitesfile: $!\n";
		while(<FILE>){
			if ($_ =~ /^>(.*)/){
				print SITES $depth.$1."\n";
			}
		}
		close SITES;
	}
	close FILE;
	
}