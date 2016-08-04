sub neva_site_entrenchment {
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..565);
	}
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[0] = $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #observed
					$hist{$bin}[1] = $total_muts/$total_length; #expected
				}
				if (!$hist{$bin}[0]){
					$hist{$bin}[0] = 0;
				}
				if (!$hist{$bin}[1]){
					$hist{$bin}[1] = 0;
				}
				print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
			}
		}
		}
	}
	
	foreach my $bin (keys %hist){
		print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	}
	
	
}

sub depth_groups_entrenchment {
	my $step = $_[0];
	my $restriction = $_[1];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;
	if ($_[2]){
		@group = @{$_[2]};
	}
	else {
		@group = (1..565);
	}
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			$self->my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_depth_hash{$ind}{$node->get_name()} > $restriction){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			#if ($total_length > 0 && $total_muts/$total_length < 0.005){
			if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					
				}
				if (!$hist{$bin}[0]){
					$hist{$bin}[0] += 0;
				}
				if (!$hist{$bin}[1]){
					$hist{$bin}[1] += 0;
				}
				#print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
			}
			}
			}
		}
		}
	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	}
	return %hist;
	
	
}

sub maxdepths_hist {
	my @array;
	my %hist;

	
	foreach my $ind (1..565){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			#print $node->get_name();
			my @args = ($ind, 1, $node);
			$self->my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			print $static_depth_hash{$ind}{$node->get_name()}."\t";
			$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			if ($total_length == 0){
				print "NA\n"
			}
			else {
				print $total_muts/$total_length."\n";
			}
		}
		}
	}
} 
 
 
 sub global_entrenchment {
	my $step = $_[0];
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..565);
	}
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;

	print "bin,observed,expected\n";
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			#print $node->get_name();
			my @args = ($ind, $step, $node);
			$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				#print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				if ($total_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
				}
			}
		}
		}
	}
	
	foreach my $bin (keys %hist){
		print $bin*$step.",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
	}
	
	
}


## full ring; each subtree has its own epsilon
sub global_entrenchment_epsilon {
	my $bin_count = $_[0];
	my @group;
	if ($_[2]){
		@group = @{$_[2]};
	}
	else {
		@group = (1..565);
	}
	my $scaled = $_[1];
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;

	print "bin,observed,expected\n";
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			print $node->get_name();
			
			my @args = ($ind, 1, $node);
			$self->my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $step = $static_depth_hash{$ind}{$node->get_name()}/$bin_count;
			print " maxdepth ".$static_depth_hash{$ind}{$node->get_name()}." step ".$step."\n";
			if ($step == 0){
				next;
			}
			@args = ($ind, $step, $node);
			$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			
			#if ($static_depth_hash{$ind}{$node->get_name()} < 100){
			my $total_muts;
			my $total_length;
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			#my $cumulative_muts;
			#my $cumulative_length;
			foreach my $bin (sort {$a <=> $b} keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				my $histbin = $bin;
				if (!$scaled){
					$histbin = $step*$bin; # real radius
				}
				
				print "bin $bin observed ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." totmut $total_muts totlen $total_length\n";
				if ($total_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#$cumulative_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
					#$cumulative_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
					#$hist{$histbin}[0] += $cumulative_muts; #observed
					#$hist{$histbin}[1] += $total_muts*$cumulative_length/$total_length; #expected
				$hist{$histbin}[0] +=  $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
				$hist{$histbin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length;
				}
			}
		}
		}
	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		print $bin.",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
	}
	
	
}


sub find_epsilon {
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..566);
	}
	my $root = $static_tree-> get_root;
	my @array;
	my @lambdas;

	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
			my @args = ($ind, 600, $node);
			$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			my $total_muts = $static_ring_hash{$ind}{$node->get_name()}{1}[0];
			if ($total_muts > 0){
				my $total_length = $static_ring_hash{$ind}{$node->get_name()}{1}[1];
				my $lambda = $total_muts;
				push @lambdas, $lambda;
				print $lambda."\n";
			}
		}
		}
	}
	
	return median(\@lambdas);

}



