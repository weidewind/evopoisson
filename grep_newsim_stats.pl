
#!/usr/bin/perl
## script for processing multiple fake results produced by FDR.sh (multiple launches of FDR_all.pl)

use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $input;
my $prot = "h1";
my $restriction = 50;
GetOptions (	
		'input=s' =>\$input,
		'prot=s' =>\$prot,
		'restriction=s' =>\$restriction,
	);

# reading simulation files
my $dirname = File::Spec->catdir(getcwd(), $input); 
my $simdirname = File::Spec->catdir($dirname, $prot); 
opendir(DH, $simdirname);
my @files = grep { /.*_[0-9]+$/ }readdir(DH); 
unless (scalar @files > 0){die "No simulation files found in folder $simdirname\n";}
closedir(DH);

my %iterations_hash;
foreach my $gulp_filename(@files){
	next if (-d $gulp_filename);
	my $fullpath = File::Spec -> catfile($simdirname, $gulp_filename);
	open GULP, "<$fullpath" or die "Cannot open $fullpath : $!";
	while (<GULP>){
				if ($_ =~ /^site/){
					# site 53 node INTNODE2434 maxdepth 269 muts 1 total_length 3512.5
					my @str_array = split(/\s+/, $_);
					$simsite = $str_array[1]; # careful 
					$simnode = $str_array[3]; # careful 
					$max_depth = $str_array[5];
					$mutnum = $str_array[7];
					$totlen = $str_array[9];
					if($max_depth > $restriction){
						push @{$iterations_hash{$simsite.":".$simnode}}, [$max_depth, $mutnum, $totlen];
					}
				}
	}
}


my $filename = $prot."_gulpselector_vector_boot_median_test_".$restriction."_single_sites";
my $filepath = File::Spec->catfile($dirname,$filename);
open FILE, "<$filepath" or die "Cannot open $filepath: $!\n";

my $sitesfile = File::Spec->catfile($dirname, $prot."_sites_with_stats");
open GROUPS, ">$sitesfile" or die "Cannot open $sitesfile: $!\n";
print GROUPS "site_node,mutations,maxlength,pvalue_epistasis(median),pvalue_epistasis(mean),pvalue_environment(median),pvalue_environment(mean),iterations,obsmedian,expmedian,obsmean,expmean,\n";
my $simsitesfile = File::Spec->catfile($dirname, $prot."_sim_sites_with_stats");
open SIMS, ">$simsitesfile" or die "Cannot open $simsitesfile: $!\n";
print SIMS "site_node,maxdepth,mutnum,totlen,bootobsmed,bootexpmed,bootobsmean,bootexpmean,\n";
print "printing to $sitesfile and $simsitesfile\n";
my $obss;
my @simstats;
while(<FILE>){
	if ($_ =~ /^[\s\t]*observed.*median/){
		my $obsmedian = (split(/:/, $_))[1]; 
		#print "obsmedian ".$obsmedian."\n";
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
	if ($_ =~ /^[\s\t]*boot[\s\t]*obs[\s\t]*median/){
		# test it
		my $bootobsmed = (split(/\s/, $_))[4];
		my $bootexpmed = (split(/\s/, $_))[8];
		$_ = <FILE>;
		$_ = <FILE>;
		my $bootobsmean = (split(/\s/, $_))[4];
		my $bootexpmean = (split(/\s/, $_))[8];
		#print $bootobsmed."\t".$bootexpmed."\t".$bootobsmean."\t".$bootexpmean."\n";
		push @simstats, [$bootobsmed, $bootexpmed, $bootobsmean, $bootexpmean];
	}

	
#	if($_ =~ /^Number of iterations/){
#		my $its = (split(/:/, $_))[1];
#		$its =~ s/[\s\t]+//g; 
#		$obss = $obss.$its;
#	}
	if($_ =~ /^No iterations found.*\s+([0-9]+)_(INTNODE[0-9]+)/){
		print SITES $1.",".$2."\n";
		$obss = "";
	}
	if ($_ =~ /^>(.*)/){
		#print $_;
		#site_node	mutations	maxlength	pvalue_epistasis(median)	pvalue_epistasis(mean)	pvalue_environment(median) iterations
		my @args = split(/[\s+]/, $1); 
		#print $1;
		foreach my $arg(@args){
			if ($arg ne ""){
				print GROUPS $arg.",";
			}
		}
		print GROUPS $obss."\n";
		$obss = "";
		foreach my $sim (@simstats){
			my $iteration_data = shift @{$iterations_hash{$args[1]}};
			print SIMS $args[1].",".$iteration_data->[0].",".$iteration_data->[1].",".$iteration_data->[2].
			",".$sim->[0].",".$sim->[1].",".$sim->[2].",".$sim->[3]."\n";
		}
		@simstats =();
	}
}
close SIMS;
close GROUPS;
close FILE;

