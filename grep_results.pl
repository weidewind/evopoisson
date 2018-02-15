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
	next unless ($filename =~ /(.*)_gulpselector_vector_boot_median_test_([0-9\.]*)_(.*)/);
	my $prot = $1;
	my $groupfile = File::Spec->catfile($dirname, $prot."_groups");
	open GROUPS, ">$groupfile" or die "Cannot open $groupfile: $!\n";
	print GROUPS "group,maxdepth,iterations,type,epi_enrichment_pvalue,env_enrichment_pvalue,epi_pvalue,env_pvalue\n";
	close GROUPS;
}	

foreach my $filename(sort @files){
	my $filepath =  File::Spec->catfile($dirname,$filename);
	next if (-d $$filepath);
	next unless ($filename =~ /(.*)_count$/);
	open FILE, "<$filepath" or die "Cannot open file $filepath: $!\n";
	my %hash;
	while(<FILE>){
		$_ =~ s/^\s++//;
		my @splitter = split(/\s++/, $_);
		if ($splitter[3]){
			$hash{$splitter[0]}{$splitter[1]} = $splitter[3];
		}
	}
	close FILE;
	my $out = $filepath."_trimmed";
	open OUT, ">$out";
	print OUT "maxdepth,group,subtree_count\n";
	foreach my $depth (sort { $a <=> $b } keys %hash){
		foreach my $group (sort keys %{$hash{$depth}}){
			print OUT $depth.",".$group.",".$hash{$depth}{$group}."\n";
			print OUT "\n" if ($group eq "all");
		}
	}
	close OUT;
}



foreach my $filename(sort @files){
	my $filepath =  File::Spec->catfile($dirname,$filename);
	next if (-d $$filepath);
	next unless ($filename =~ /(.*)_gulpselector_vector_boot_median_test_([0-9\.]*)_(.*)/);
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
				my @splitter = split(/[\s\t]+/, $_);
				print GROUPS $group.",".$depth.",".$its.",";
				if (scalar @splitter == 3){
					print GROUPS $splitter[0].",,,".$splitter[1].",".$splitter[2]."\n";
				}
				elsif(scalar @splitter == 5){
					print GROUPS $splitter[0].",".$splitter[1].",".$splitter[2].",".$splitter[3].",".$splitter[4]."\n";
				}
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
						my @args = split(/[:\s+]/, $1); 
						pop @args; # iterations number
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
	else {
		$filename =~ /(.*)_gulpselector_vector_boot_median_test_([0-9\.]*)_(.*)/;
		my $group = $3;
		my $maxdepth = $2;
		my $groupfile = File::Spec->catfile($dirname, $prot."_shift_groups");
		open GROUPS, ">>$groupfile" or die "Cannot open $groupfile: $!\n";
		my $obsmean;
		my $expmean;
		my $obsmedian;
		my $expmedian;
		while(<FILE>){
			if ($_ =~ /^\s+observed\s+mean/){
				$obsmean = (split(/:/, $_))[1];
				$obsmean =~ s/[\s\t]+//g; 
			}
			if ($_ =~ /^\s+observed\s+median/){
				$obsmedian = (split(/:/, $_))[1];
				$obsmedian =~ s/[\s\t]+//g; 
			}
			if ($_ =~ /^\s+poisson\s+expected\s+mean/){
				$expmean = (split(/:/, $_))[1];
				$expmean =~ s/[\s\t]+//g; 
			}
			if ($_ =~ /^\s+poisson\s+expected\s+median/){
				$expmedian = (split(/:/, $_))[1];
				$expmedian =~ s/[\s\t]+//g; 
			}
		}
		my $meandiff = $obsmean-$expmean; 
		my $mediandiff = $obsmedian-$expmedian; 
		print GROUPS $group.",".$maxdepth.",".$meandiff.",".$mediandiff."\n";
		close GROUPS;
	}
	close FILE;
	
}
