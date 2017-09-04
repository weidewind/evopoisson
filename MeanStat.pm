package MeanStat;
use base ("AgeingStat");

sub new {
    my $class = shift;
    return bless {}, $class;
}

sub computeStats{
	return $self->computeDiff($args);
}


sub hist_stat {
	my $self = shift;
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


######################## not used any more

#takes a hash of probabilities for any distances (including nonintegral)
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



sub hist_mean_for_hash_arr{
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

	return hist_mean(\@hist);
}

