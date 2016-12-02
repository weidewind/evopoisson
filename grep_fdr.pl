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
			next unless ($filename =~ /(.*)_gulpselector_vector_boot_median_test_([0-9]+)_(.*)/);
			my $filepath = File::Spec->catfile($dipath,$filename);
			my $prot = $1;
			my $depth = $2;
			my $ending = $3;
			open FILE, "<$filepath" or die "Cannot open $filename: $!\n";
			if ($ending =~ /(all)$/){
				my $group = $1;
				my $groupfile = File::Spec->catfile($dirname, $prot."_allnew_fdr");
				open GROUPS, ">>$groupfile" or die "Cannot open $groupfile: $!\n";
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
						$obss = $obss.$expmean;
					}
					if ($_ =~ /^me/){
						chomp($_);
						$_ =~ s/[\s\t]+/,/g;
						print GROUPS $di.",".$group.",".$depth.",".$_.",".$obss."\n";
						my ($stat, $epi, $env) = split(/\s+/, $_);
						#print "stat ".$stat." epi ".$epi." env ".$env."\n";
						push @{$hash{$prot}{$depth}{$ending}{$stat}{"epi"}}, $epi;
						push @{$hash{$prot}{$depth}{$ending}{$stat}{"env"}}, $env;
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
	}
}
	my @thresholds = (0.05, 0.01, 0.005, 0.001, 0.0005);
	my @stats = ("median_stat", "mean_stat");
	foreach my $prot (keys %hash){
		foreach my $depth (keys %{$hash{$prot}}){
			foreach my $ending (keys %{$hash{$prot}{$depth}}){
				foreach my $stat (keys %{$hash{$prot}{$depth}{$ending}}){
					foreach my $thr (@thresholds){
						my $counter;
						my $totals;
						foreach my $dot(@{$hash{$prot}{$depth}{$ending}{$stat}{"epi"}}){
							#print "dot ".$dot."\n";
							$totals++;
							if ($dot <= $thr){$counter++;}
						}
						print $prot." ".$depth." ".$ending." ".$stat." ".$thr." epi ".$counter/$totals."\n";
						my $counter;
						my $totals;
						foreach my $dot(@{$hash{$prot}{$depth}{$ending}{$stat}{"env"}}){
							#print "dot ".$dot."\n";
							$totals++;
							if ($dot <= $thr){$counter++;}
						}
						print $prot." ".$depth." ".$ending." ".$stat." ".$thr." env ".$counter/$totals."\n";
						
						}
				}
			}
		}
	}
