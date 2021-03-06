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
		my $its;
		while(<FILE>){
			if($_ =~ /^Number of iterations/){
						$its = (split(/:/, $_))[1];
						$its =~ s/[\s\t]+//g; 
			}
			if ($_ =~ /^me/){
				print GROUPS $group."\t".$depth."\t".$its."\t".$_;
			}
			if ($_ =~ /^#/){
				$its = "";
			}
		}
		close GROUPS;
		 
	}
	elsif ($ending =~ /.*single_sites$/){
		my $sitesfile = File::Spec->catfile($dirname, $prot."_sites");
		open SITES, ">>$sitesfile" or die "Cannot open $sitesfile: $!\n";
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
					if($_ =~ /^Number of iterations/){
						my $its = (split(/:/, $_))[1];
						$its =~ s/[\s\t]+//g; 
						$obss = $obss.$its;
					}
					if($_ =~ /^No iterations found.*\s+([0-9]+)_(INTNODE[0-9]+)/){
						print SITES $depth.",".$1.",".$2."\n";
						$obss = "";
					}
					if ($_ =~ /^>(.*)/){
						my @args = split(/[_\s+]/, $1); 
						print SITES $depth.",";
						foreach my $arg(@args){
							if ($arg ne ""){
								print SITES $arg.",";
							}
						}
						print SITES $obss."\n";
						$obss = "";
						
					}

				
		}
		close SITES;
	}
	close FILE;
	
}