#!/usr/bin/perl
package Enrichment;
require Exporter;	
@ISA = qw(Exporter);	
@EXPORT = qw(get_groups only_groups parse_real_stats diffvalues get_fake_dirs parse_fake_stats meanvalues);
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Mutmapnolim;
use Groups;
use Cwd qw(abs_path cwd getcwd);
use Sort::Rank  qw(rank_sort);
use Data::Dumper;
use List::Util  qw(sum min max);


sub diffvalues{
	my @sites = @{$_[0]};
	my $groups = $_[1];
	my $group_names = $_[2];
	my $value = $_[3];
	my ($g, $gn) = only_groups($groups, $group_names);
	my ($c, $cn) = only_compl($groups, $group_names);
	my %realgranks = meanvalues(\@sites,$g, $gn, $value);
	my %realcranks = meanvalues(\@sites,$c, $cn, $value);
	my %diff;
	foreach my $name (keys %realgranks){
		$diff{$name} = $realgranks{$name}-$realcranks{$name};
	}	
	return %diff;
}
sub get_groups{
	my $data = shift;
	my $protein = shift;
	$fastafile = File::Spec->catfile(getcwd(), $data, $protein.".all.fa");
	my @arr = Mutmapnolim::parse_fasta($fastafile);
	my $alignment_length = $arr[1];
	my $length = $alignment_length/3;
	my @groups_and_names = Groups::get_predefined_groups_and_names_for_protein($protein, $length);
	my @groups = @{$groups_and_names[0]};
	my @group_names = @{$groups_and_names[1]};
	return \@groups, \@group_names;
}




sub parse_real_stats{
	my $input = shift;
	my $protein = shift;
	my $dir = File::Spec->catdir(getcwd(), $input);
	my $sitesfile = File::Spec->catfile($dir, $protein."_sites");
	open SITES, $sitesfile or die "Cannot open $sitesfile: $!";
	my @sites;
	while(<SITES>){
		my @splitter = split(/,/);
		my $row;
		if (scalar @splitter == 16){
			$row = {     site =>	$splitter[1],
						 node =>	$splitter[2],
						 mutations  =>	$splitter[3], 	
						 maxlength =>	$splitter[4],
						 'pvalue_epistasis(mean)' =>	$splitter[5],
						 'pvalue_epistasis(median)' => $splitter[6],
						 'pvalue_epistasis(bp)' => $splitter[7],
						 'pvalue_environment(mean)' => $splitter[8],
						 'pvalue_environment(median)' => $splitter[9],	
						 'pvalue_environment(bp)' => $splitter[10],
						 obsmean => $splitter[11],
						 expmean => $splitter[12],
						 obsmedian => $splitter[13],
						 expmedian => $splitter[14],
						 iterations => $splitter[15]};
		}
		else {die "Incorrect sites file format!";}
		push @sites, $row;
	}
	close SITES;
	return @sites;
}

sub parse_fake_stats{
	my $input = shift;
	my $protein = shift;
	my $dir = File::Spec->catdir(getcwd(), $input);
	my $sitesfile = File::Spec->catfile($dir, $protein."_sites_with_stats");
	open SITES, $sitesfile or die "Cannot open $sitesfile: $!";
	my @sites;
	my $h = <SITES>;
	while(<SITES>){
		my @splitter = split(/[,:]/);
		my $row = {      site =>	$splitter[0],
						 node =>	$splitter[1],
						 mutations  =>	$splitter[2], 	
						 maxlength =>	$splitter[3],
						 'pvalue_epistasis(mean)' =>	$splitter[4],
						 'pvalue_epistasis(median)' => $splitter[5],
						 'pvalue_epistasis(bp)' => $splitter[6],
						 'pvalue_environment(mean)' => $splitter[7],
						 'pvalue_environment(median)' => $splitter[8],	
						 'pvalue_environment(bp)' => $splitter[9],
						 iterations => $splitter[10]
		};
		push @sites, $row;
	}
	close SITES;
	return @sites;
}

sub get_fake_dirs {
	my $fakesdir = shift;
	my $tail = shift;
	my $dirname = File::Spec->catdir(getcwd(), $fakesdir); 
	opendir(DH, $dirname);
	my @dirs = readdir(DH);
	closedir(DH);
	if (scalar @dirs == 0){
		die "No files or dirs found in $input!\n";
	}
	#print "dirname ".$dirname."\n";
	my @fulldirs;
	foreach my $di(sort @dirs){
		my $dipath =  File::Spec->catdir($dirname,$di);
	#	print "dipath ".$dipath."\n";
		if (-d $dipath && $di =~ /([0-9]+)_fake/){
			my $inside =  File::Spec->catdir($fakesdir,$di,$tail);
			push @fulldirs, $inside;
		#	print "inside ".$inside."\n";
		}
	}	
	return @fulldirs;
}



sub meanranks{
	my @sites = @{$_[0]};
	my @groups = @{$_[1]};
	my @group_names = @{$_[2]};
	my $value = $_[3];
	my $action = $_[4];
	
	my @sorted = rank_sort(\@sites, sub {
	        # Extract score from an element
	        my $item = shift;
	        return 1-$item->{$value};
	    });
	
	my %groupranks;
	for (my $i = 0; $i <scalar @groups; $i++){
		my $groupname = $group_names[$i]; 
	#	print $groupname."\n";
		my @ranks;
		foreach my $s(@{$groups[$i]}){
			foreach my $elm(@sorted){
				if ($elm->[2]->{site} eq $s){
					push @ranks, $elm->[0];
				}
			}
		}
		if ($action) {
			$groupranks{$groupname} = &$action(\@ranks);
		}
		else {$groupranks{$groupname} = sum(@ranks)/(scalar @ranks)};
	}   
	
#	print Dumper(\%groupranks);
	return %groupranks;
}

sub meanvalues{
	my @sites = @{$_[0]};
	my @groups = @{$_[1]};
	my @group_names = @{$_[2]};
	my $value = $_[3];
	my $action = $_[4];
	
	my %groupranks;
	for (my $i = 0; $i <scalar @groups; $i++){
		my $groupname = $group_names[$i]; 
		my @ranks;
		foreach my $s(@{$groups[$i]}){
			foreach my $elm(@sites){
				if ($elm->{site} eq $s){
					push @ranks, $elm->{$value};
				}
			}
		}
		if (! @ranks) {$groupranks{$groupname} = undef;}
		else {
			if ($action) {
				$groupranks{$groupname} = &$action(\@ranks);
			}
			else {$groupranks{$groupname} = sum(@ranks)/(scalar @ranks)};
		}
	}   
	return %groupranks;
}

sub only_groups {
	my @groups = @{$_[0]};
	my @group_names = @{$_[1]};
	my @g;
	my @n;
	for (my $i = 0; $i <scalar @groups; $i= $i+2){
		push @g, $groups[$i];
		push @n, $group_names[$i];
	}
	return (\@g, \@n);
}

sub only_compl {
	my @groups = @{$_[0]};
	my @group_names = @{$_[1]};
	my @g;
	my @n;
	for (my $i = 1; $i <scalar @groups; $i= $i+2){
		push @g, $groups[$i];
		push @n, $group_names[$i-1];
	}
	push @g, [];
	push @n, "all";
	return (\@g, \@n);
}
1;
