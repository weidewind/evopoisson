#!/usr/bin/perl
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $input;
my $state = "nsyn";
my $no_neighbour_changing;
GetOptions (	
		'input=s' =>\$input,
		'state=s' =>\$state,
		'no_neighbour_changing' =>\$no_neighbour_changing,
	);

my $dirname = File::Spec->catdir(getcwd(), $input); 

opendir(DH, $dirname);
my @dirs = readdir(DH);
closedir(DH);
if (scalar @dirs == 0){
	die "No files or dirs found in $input!\n";
}
my $outfile = File::Spec->catfile($dirname, "fake_fdr_stats");
open GROUPS, ">$outfile" or die "Cannot open $outfile: $!\n";
print GROUPS "fakeno,group,maxdepth,subtree_count,iterations,type,epi_enrichment_pvalue,env_enrichment_pvalue,epi_pvalue,env_pvalue,stat\n";
foreach my $di(sort @dirs){
	my $dipath =  File::Spec->catdir($dirname,$di);
	if (-d $dipath && $di =~ /([0-9]+)_fake/){
		my $fakeno = $1;
		
		$dipath =  File::Spec->catdir($dipath,$state,"maxpath_not_subtracted");
		if ($no_neighbour_changing){
			$dipath =  File::Spec->catdir($dipath,"no_neighbour_changing");
		}
		opendir(DH, $dipath);
		print "fake dir found: $dipath \n";
		my @files = readdir(DH);
		closedir(DH);
		if (scalar @files == 0){
			warn "No files found in $dipath!\n";
			next;
		}
		
		
		my ($countfile) = grep(/(.*)_count/, @files);	
		my %hash;
		next unless defined $countfile;
		my $countpath =  File::Spec->catfile($dipath,$countfile);
		open COUNT, "<$countpath" or die "Cannot open $countpath: $!\n";
		while (<COUNT>){
			my @splitter = split (/[\s\t]+/, $_);
			if ($_ =~ /^T/){
				$hash{$splitter[3]}{$splitter[2]} = $splitter[4]; # hash{group_name}{restriction} = number of subtrees
			}
			else {
				$hash{$splitter[1]}{$splitter[0]} = $splitter[3];
			}
			
		}
		close COUNT;
		foreach my $filename(sort @files){
			#print "$filename found\n";
			my $filepath =  File::Spec->catfile($dipath,$filename);
			next if (-d $filepath);
			next unless ($filename =~ /(.*)_gulpselector_vector_boot_median_test_([0-9\.]+)_(.*)/);
			my $prot = $1;
			my $depth = $2;
			my $ending = $3;
			if ($ending =~ /(all)$/){
				print "Going to open $filepath \n";
			}

			my $line;
			open FILE, "<$filepath" or die "Cannot open $filepath: $!\n";
			if ($ending =~ /(.*)_complement$/ || $ending =~ /(all)$/){
				my $group = $1;
				my $its;
				my $median_stat;
				my $mean_stat;
				my $median_line;
				my $mean_line;
				while(<FILE>){
					if($_ =~ /^Number of iterations/){
								$its = (split(/:/, $_))[1];
								$its =~ s/[\s\t]+//g; 
					}
					if ($ending =~ /(all)$/){
						if ($_ =~ /^[\s\t]+observed\smean/){
							my @splitter = split(/[\s\t]+/, $_);
							my $obsmean = pop @splitter;
							$mean_stat += $obsmean;
						}
						if ($_ =~ /^[\s\t]+observed\smedian/){
							my @splitter = split(/[\s\t]+/, $_);
							my $obsmed = pop @splitter;
							$median_stat += $obsmed;
						}
						if ($_ =~ /^[\s\t]+poisson\sexpected\smean/){
							my @splitter = split(/[\s\t]+/, $_);
							my $expmean = pop @splitter;
							$mean_stat -= $expmean;
						}
						if ($_ =~ /^[\s\t]+poisson\sexpected\smedian/){
							my @splitter = split(/[\s\t]+/, $_);
							my $expmedian = pop @splitter;
							$median_stat -= $expmedian;
						}
					}

					if ($_ =~ /^me/){
						my @splitter = split(/[\s\t]+/, $_);
						my $type = shift @splitter;
						my $line = $type;
						if ($group eq "all"){
							$line = $line.",,,".$splitter[0].",".$splitter[1];
						}
						else {
							$line = $line.",".$splitter[0].",".$splitter[1].",".$splitter[2].",".$splitter[3];
						}
						if ($ending =~ /(all)$/){
							if ($_ =~ /^median/){$line = $line.",".$median_stat;}
							if ($_ =~ /^mean/){$line = $line.",".$mean_stat;}
							$line = $line."\n";
							print GROUPS $fakeno.",".$group.",".$depth.",".$hash{$group}{$depth}.",".$its.",".$line;
							
						}
						else {
							if ($type eq "median_stat"){
								$median_line = $fakeno.",".$group.",".$depth.",".$hash{$group}{$depth}.",".$its.",".$line;
							}
							if ($type eq "mean_stat"){
								$mean_line = $fakeno.",".$group.",".$depth.",".$hash{$group}{$depth}.",".$its.",".$line;
							}
						}
						#print GROUPS $fakeno.",".$group.",".$depth.",".$hash{$group}{$depth}.",".$its.",".$line;
					}
					if ($_ =~ /^#/){
						$its = "";
					}
				}
				close FILE;
				if ($ending =~ /(.*)_complement$/){
					my $groupfilename = $prot."_gulpselector_vector_boot_median_test_".$depth."_".$group;
					my $groupfilepath =  File::Spec->catfile($dipath,$groupfilename);
					open FILE, "<$groupfilepath" or die "Cannot open $groupfilepath: $!\n";
					while (<FILE>){
						if ($_ =~ /^[\s\t]+observed/){
							my @splitter = split(/[\s\t]+/, $_);
							print $splitter[3]." ".$splitter[7]." ".$splitter[10]." ".$splitter[14]."\n";
							$median_stat = $splitter[3]-$splitter[7];
							$mean_stat = $splitter[10]-$splitter[14];
							last;
						}
					}
					close FILE;
					print GROUPS $median_line.",".$median_stat."\n";
					print GROUPS $mean_line.",".$mean_stat."\n";
				}
				 	
			}
		
			
		
			
			
		}
	}
}	
close GROUPS;	