package MedianStat;
use base ("AgeingStat");

sub new {
    my $class = shift;
    return bless {}, $class;
}

sub computeStats{
	return $self->computeDiff($args);
}


#takes an array of probabilities for 0,1,2...
sub hist_stat{
	my $self = shift;
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
		if (abs ($head - $summ/2) < 0.000001) {last;}
	}
	if (abs ($head - $summ/2) < 0.000001){ # magic 
		my $leftmedian = $median;
		my $rightmedian;
		my $newhead = $head;
		while ($newhead == $head){
			$newhead += $hist[$interval];
			$rightmedian = $interval*$step;
			$interval++;
		}
		$median = ($leftmedian+$rightmedian)/2;
	}
	return $median;
}




############ not used anymore


#takes a hash of probabilities for any distances (including nonintegral)
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

sub hist_median_for_hash_arr{
	my %prehist =  %{$_[0]};
	my $number = $_[1];
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i}[$number];
		}
		else {
			push @hist, 0;
		}
	}
	return hist_median(\@hist);
}