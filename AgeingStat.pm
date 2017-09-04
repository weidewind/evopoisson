package AgeingStat;
use Switch;
use List::Util qw(sum min max);

sub new {
	my ($class, $type) = @_;
	switch ($type) {
		case "mean"		{ $class = "MeanStat"; }
		case "median"	{ $class = "MedianStat";  }
		case "bp"	{ $class = "BPStat";  }
	}
	return bless {}, $class;			
}

sub computeDiff {
	my $self = shift;
	my ($args) = @_;
	my $obs = $self->hist_stat_for_hash($args->{obshash},$args->{step});
	my $exp = $self->hist_stat_for_hash($args->{exphash},$args->{step});
	my $value = $obs-$exp;
	$self->{'obs'} = $obs;
	$self->{'exp'} = $exp;
	$self->{'value'} = $value;
	return $value;
}

#takes a hash of probabilities for 0,1,2...
sub hist_stat_for_hash{
	my $self = shift;
	my @hist = hist_to_array($_[0]);
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	if (! defined $step) {$step  = 1;}
	return $self->hist_stat(\@hist, $step);
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

1;