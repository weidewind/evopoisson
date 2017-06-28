#!/usr/bin/perl
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $input;
my $state = "nsyn";
my $no_neighbour_changing;
my $skip_stoppers;
GetOptions (	
		'input=s' =>\$input,
		'state=s' =>\$state,
		'no_neighbour_changing' =>\$no_neighbour_changing,
		'skip_stoppers' =>\$skip_stoppers,
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
		$dipath =  File::Spec->catdir($dipath,$state,"maxpath_not_subtracted");
		if ($no_neighbour_changing){
			$dipath =  File::Spec->catdir($dipath,"no_neighbour_changing");
		}
		if($skip_stoppers){
			$dipath =  File::Spec->catdir($dipath,"skip_stoppers");
		}
		print "fake dir found: $dipath \n";
		opendir(DH, $dipath);
		my @files = readdir(DH);
		closedir(DH);
		if (scalar @files == 0){
			warn "No files found in $dipath!\n";
			next;
		}
		
		foreach my $filename(sort @files){

			my $filepath =  File::Spec->catfile($dipath,$filename);
			next if (-d $filepath);
			if ($filename =~ /(.*)_sim_sites_with_stats/){
				print "$filename found\n";
				open FILE, "<$filepath" or die "Cannot open $filepath: $!\n"; 
				my $str =  <FILE>; 
				while(<FILE>){
					print SIMS $fakeno.",".$_;
				}
				close FILE;
				
			}
			if ($filename =~ /[a-z][0-9]_sites_with_stats/){
				print "$filename found\n";
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
		
		
		
		
		
		
	