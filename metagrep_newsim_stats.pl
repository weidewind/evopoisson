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
my @dirs = readdir(DH);
closedir(DH);
if (scalar @dirs == 0){
	die "No files or dirs found in $input!\n";
}
my $sitesfile = File::Spec->catfile($dirname, "meta_sites_with_stats");
open SITES, ">$sitesfile" or die "Cannot open $sitesfile: $!\n";
print SITES "fakeno,site_node,mutations,maxlength,pvalue_epistasis(median),pvalue_epistasis(mean),pvalue_environment(median),pvalue_environment(mean),iterations,obsmedian,expmedian,obsmean,expmean,\n";

my $simfile = File::Spec->catfile($dirname, "meta_sim_sites_with_stats");
open SIMS, ">$simfile" or die "Cannot open $simfile: $!\n";
print SIMS "fakeno,site_node,maxdepth,mutnum,totlen,bootobsmed,bootexpmed,bootobsmean,bootexpmean,\n";
foreach my $di(sort @dirs){
	my $dipath =  File::Spec->catdir($dirname,$di);
	if (-d $dipath && $di =~ /([0-9]+)_fake/){
		my $fakeno = $1;
		print "fake dir found: $di \n";
		$dipath =  File::Spec->catdir($dipath,"nsyn","maxpath_not_subtracted");
		opendir(DH, $dipath);
		my @files = readdir(DH);
		closedir(DH);
		if (scalar @files == 0){
			die "No files found in $input!\n";
		}
		
		foreach my $filename(sort @files){
			print "$filename found\n";
			my $filepath =  File::Spec->catfile($dipath,$filename);
			next if (-d $filepath);
			if ($filename =~ /(.*)_sim_sites_with_stats/){
				open FILE, "<$filepath" or die "Cannot open $filepath: $!\n"; 
				my $str =  <FILE>; 
				while(<FILE>){
					print SIMS $fakeno.",".$_;
				}
				close FILE;
				
			}
			if ($filename =~ /[a-z][0-9]_sites_with_stats/){
				open FILE, "<$filepath" or die "Cannot open $filepath: $!\n";  
				my $str =  <FILE>; 
				while(<FILE>){
					print SITES $fakeno.",".$_;
				}
				close FILE;
			}
		}
	}
}			
		
		
		
		
		
		
	