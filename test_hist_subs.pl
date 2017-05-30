#!/usr/bin/perl

#use Mutmapnolim qw(hist_median_for_hash hist_mean_for_hash hist_median hist_mean hist_to_array);
use List::Util qw(sum);
use strict;


sub hist_median_for_hash{
	my @hist = hist_to_array($_[0]);
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	my $verbose = $_[2];
	if (! defined $step) {$step  = 1;}
	return hist_median(\@hist, $step, $verbose);
}


#takes a hash of probabilities for 0,1,2...
sub hist_mean_for_hash {
	my @hist = hist_to_array($_[0]);
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	if (! defined $step) {$step  = 1;}
	return hist_mean(\@hist, $step);
}



sub hist_median_for_hash_nobin{
    my %hist = %{$_[0]};
	my $verbose = $_[2];
	my $summ = sum (values %hist);
	my $head = 0;
	my @distances = sort { $a <=> $b } keys %hist;
	my $distance = shift @distances;
	my $median = 0;
	
	while ($head < $summ/2){
		unless ($hist{$distance} == 0) {$head += $hist{$distance};}
		$median = $distance; 
		$distance = shift @distances;
		if ($verbose) {
			my $half = $summ/2;
			print $head." ".$median." ".$half."\n";
		}
		if (abs ($head - $summ/2) < 0.000001) {last;}
	}
	if ($verbose) {print $summ."\n";}
	if (abs ($head - $summ/2) < 0.000001){ # magic 
	#	$median += 0.5*$step;
		my $leftmedian = $median;
		my $rightmedian;
		my $newhead = $head;
		if ($verbose) {print "head == summ/2: ".$head." ".$leftmedian."\n";}
		while ($newhead == $head){
			$newhead += $hist{$distance};
			$rightmedian = $distance;
			if ($verbose) {print $newhead." ".$rightmedian."\n";} 
			$distance = shift @distances;
		}
		$median = ($leftmedian+$rightmedian)/2;
		if ($verbose) {print "median is $median\n";}
	}
#print_hist(\@hist);
	if ($verbose) {print "median is $median\n";}
	return $median;
}

#takes a hash of probabilities for 0,1,2...
sub hist_mean_for_hash_nobin {
	my %hist = %{$_[0]};
	my $summ =  sum (values %hist);
	my @distances = sort { $a <=> $b } keys %hist;
	my $integer;
	for my $dist(@distances){
		$integer += $dist*$hist{$dist};
	}
	if ($summ > 0){
		return $integer/$summ;
	}
	else {
		return "NaN";
	}
}


#takes an array of probabilities for 0,1,2...
sub hist_median{
	my @hist = @{$_[0]};
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	my $verbose = $_[2];
	if (! defined $step) {$step  = 1;}
	my $summ = sum (@hist);
	my $head = 0;
	my $interval = 0;
	my $median = 0;
	
	while ($head < $summ/2){
		unless ($hist[$interval] == 0) {$head += $hist[$interval];}
		$median = $interval*$step; 
		$interval++;
	#	if ($verbose) {
	#		my $half = $summ/2;
	#		print $head." ".$median." ".$half."\n";
	#	}
		if (abs ($head - $summ/2) < 0.000001) {last;}
	}
	#if ($verbose) {print $summ."\n";}
	if (abs ($head - $summ/2) < 0.000001){ # magic 
	#	$median += 0.5*$step;
		my $leftmedian = $median;
		my $rightmedian;
		my $newhead = $head;
	#	if ($verbose) {print "head == summ/2: ".$head." ".$leftmedian."\n";}
		while ($newhead == $head){
			$newhead += $hist[$interval];
			$rightmedian = $interval*$step;
		#	if ($verbose) {print $newhead." ".$rightmedian."\n";} 
			$interval++;
		}
		$median = ($leftmedian+$rightmedian)/2;
	#	if ($verbose) {print "median is $median\n";}
	}
#print_hist(\@hist);
#	if ($verbose) {print "median is $median\n";}
	return $median;
}

sub hist_mean {
	my @hist = @{$_[0]};
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	if (! defined $step) {$step  = 1;}
	my $summ = sum (@hist);
	my $integer;
	for(my $i = 0; $i <scalar @hist; $i++){
		$integer += $i*$hist[$i]*$step;
	}
	if ($summ > 0){
		return $integer/$summ;
	}
	else {
		return "NaN";
	}
}



# takes hash of probabilities for 0,1,2... and returns an array of ordered values
sub hist_to_array {
	my %prehist =  %{$_[0]};
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	# for (my $i = 1; $i <= $sorted_keys[-1]; $i++){ # changed at 21.09.2016
	for (my $i = 0; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i};
		}
		else {
			push @hist, 0;
		}
	}
	return @hist;
}


my %hist = (0 => 2, 1 => 2, 2 => 2, 3 => 2);
test(\%hist);

my %hist = (0 => 3, 1 => 1, 2 => 1, 3 => 1);
test(\%hist);

my %hist = (0 => 0.5, 1 => 1, 2 => 1.5, 3 => 4);
test(\%hist);

my %hist = (0 => 0.5, 1 => 1, 2 => 1.5, 3 => 4, 4 => 6);
test(\%hist);

my %hist = (0 => 2.5, 1 => 1, 2 => 1.5, 3 => 4, 4 => 1, 10 => 2);
test(\%hist);

my %hist = (22 => 8, 0 => 2.5, 1 => 1, 2 => 1.5, 3 => 4, 4 => 1, 10 => 2);
test(\%hist);

my %hist = (0.5 => 0.5, 2.7 => 1.5, 3 => 4, 10 => 4);
test(\%hist);
my @sorted_keys = sort keys %hist;
print @sorted_keys;
sub test {
	my %hist = %{$_[0]};
	my $prev_med = hist_median_for_hash(\%hist);
	my $prev_mean  = hist_mean_for_hash(\%hist);
	print $prev_med."\t".$prev_mean."\n";
	
	my $new_med = hist_median_for_hash_nobin(\%hist);
	my $new_mean  = hist_mean_for_hash_nobin(\%hist);
	print $new_med."\t".$new_mean."\n";
	
	if ($new_med != $prev_med ||$new_mean != $prev_mean){print "Error!\n"}
}




