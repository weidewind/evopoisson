package BPStat;
use base ("AgeingStat");
use List::Util qw(sum min max);
use Data::Dumper;

sub new {
    my $class = shift;
        my $self = {
    	name => "bp",
    };
    return bless $self, $class;	
}

sub printStats {
	my $self = shift;
	my $output = shift;
	my $str = "\n".$self->{name}." value: ".$self->{value}."\n";
	if ($output) { print $output $str;}
	else {print $str;}
}

sub computeStats{
	my $self = shift;
	my ($args) = @_;
	my @plot = tttplot($args->{obshash}, $args->{exphash});
	my $value = BPstat(\@plot, $args->{zscore});
	$self->{'value'} = $value;
	return $value;
}

sub tttplot {
	my $obshist = $_[0];
	my $exphist = $_[1];
	
   	my @allbins = keys %{{map {($_ => 1)} (keys %{$exphist}, keys %{$obshist})}};
	my @sortedbins = sort {$a <=> $b} @allbins;
	my $totmuts = sum(values %{$obshist});
	
	my @ttt;
	my $ti;
	my $ni;
	for my $bin (@sortedbins){
		$ti += $exphist->{$bin};
		for (my $muts = 0; $muts < $obshist->{$bin}; $muts++){
			$ni += 1;
			push @ttt, [$ni/$totmuts, $ti/$totmuts];
		}
	}
	push @ttt, [$ni/$totmuts, $ti/$totmuts];
	return @ttt;
}

sub norm_plot{
	my @plot = @{$_[0]};
	foreach my  $dot(@plot){
		$dot->[1] = $dot->[1]/$plot[-1][1];
	} 
	return @plot;
}

sub BPstat {
	my @tttplot = @{$_[0]};
	my $zscore = $_[1];
	my $lastdot;
	if ($tttplot[0][0] == 1){
			print "Single mutation, cannot compute BP statistics!\n";
			return undef;
	}
	print "---------\n";
	foreach my $dot(@tttplot){
		print $dot->[0].",".$dot->[1]."\n";
	}
	while ($tttplot[-1][0] == 1){$lastdot = pop @tttplot;}
	my $w;
	foreach my $dot(@tttplot){
		$dot->[1] = $dot->[1]/$lastdot->[1];
		$w += $dot->[1];
	} 
	if ($zscore){
		my $nminus1 = scalar @tttplot;
		my $z = ($w - $nminus1/2)/($nminus1/12)**0.5;
		return $z;
	}
	else {
		return $w;
	}
}


1;