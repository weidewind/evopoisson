#!/usr/bin/perl
## script for processing multiple fake results produced by FDR.sh (multiple launches of FDR_all.pl)

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

my %hash;
my $sitesfile = File::Spec->catfile($dirname, "h1_sites_fdr_with_stats");
open GROUPS, ">>$sitesfile" or die "Cannot open $sitesfile: $!\n";
print GROUPS "depth,site,node,mutations,maxlength,pvalue_epistasis_median,pvalue_epistasis_mean,pvalue_environment_median,pvalue_environment_mean,iterations,obs_median,exp_median,obs_mean,exp_mean\n";
close GROUPS;
foreach my $di(sort @dirs){
	
	my $dipath =  File::Spec->catdir($dirname,$di);
	if (-d $dipath && $di =~ /_fake_(.*)/){
		#my $num = $1;
		#unless ($num <= 460) {next;} # not all folders contain new results
		print "fake dir found: $di \n";
		opendir(DI, $dipath);
		my @filenames = readdir(DI);
		close DI;
		foreach my $filename(sort @filenames){
			next unless ($filename =~ /(.*)_gulpselector_vector_boot_median_test_([0-9]+)_single_sites/);
			my $filepath = File::Spec->catfile($dipath,$filename);
			my $prot = $1;
			my $depth = $2;
			open FILE, "<$filepath" or die "Cannot open $filename: $!\n";

				my $sitesfile = File::Spec->catfile($dirname, $prot."_sites_fdr_with_stats");
				open GROUPS, ">>$sitesfile" or die "Cannot open $sitesfile: $!\n";
				my $obss;
				while(<FILE>){
					if ($_ =~ /^[\s\t]*observed.*median/){
						my $obsmedian = (split(/:/, $_))[1]; 
						#print $obsmedian."\n";
						$obsmedian =~ s/[\s\t]+//g;
						#print $obsmedian."\n";
						$obss = $obss.$obsmedian.",";
					}
					if ($_ =~ /^[\s\t]*poisson.*median/){
						my $expmedian = (split(/:/, $_))[1]; 
						$expmedian =~ s/[\s\t]+//g;
						$obss = $obss.$expmedian.",";
					}
					if ($_ =~ /^[\s\t]*observed.*mean/){
						my $obsmean = (split(/:/, $_))[1]; 
						$obsmean =~ s/[\s\t]+//g; 
						$obss = $obss.$obsmean.",";
					}
					if ($_ =~ /^[\s\t]*poisson.*mean/){
						my $expmean = (split(/:/, $_))[1];
						$expmean =~ s/[\s\t]+//g; 
						$obss = $obss.$expmean.",";
					}
				#	if($_ =~ /^Number of iterations/){
				#		my $its = (split(/:/, $_))[1];
				#		$its =~ s/[\s\t]+//g; 
				#		$obss = $obss.$its;
				#	}
					if($_ =~ /^No iterations found.*\s+([0-9]+)_(INTNODE[0-9]+)/){
						print SITES $depth.",".$1.",".$2."\n";
						$obss = "";
					}
					if ($_ =~ /^>(.*)/){
						my @args = split(/[_\s+]/, $1); 
						print GROUPS $depth.",";
						foreach my $arg(@args){
							if ($arg ne ""){
								print GROUPS $arg.",";
							}
						}
						print GROUPS $obss."\n";
						$obss = "";
						
					}

				}
				close GROUPS;
				 
			

			close FILE;
		}
	}
}