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

## used up to 4.02. site_node is node_name only, wrong FILE output
sub depth_groups_entrenchment_optimized_selection_alldepths_before402 {
	my $step = $_[0];
	my $ancestral_nodes = $_[2];
	my $subtract_tallest = $_[3];
	my $tag = $_[4];
	
	my $filename = $static_protein."_for_enrichment_".$tag;
	open FILE, ">>$filename";
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print FILE ">new iteration\n";
	
	my @group;
	if ($_[5]){
		@group = @{$_[5]};
	}
	else {
		@group = (1..565);
	}
	
	
	foreach my $key (keys %{$ancestral_nodes}){
	#print " ancestral node $key \n";
	}
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root, $subtract_tallest);
	$self->my_visit_depth_first ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);

	foreach my $ind (@group){
	#print "here my ind $ind\n";
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $site_node = $node->get_name();
			#print "here my site_node $site_node\n";
			if (!$ancestral_nodes->{$site_node}){
				#print "$ind $site_node NOT IN REAL ANCESTORS\n";
				my $totmut;
				my $totlen;
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$totmut += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$totlen += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
				}	
			#	print " totmut $totmut, totlen $totlen\n";
				next;
			}
		#	print "$site_node is still here\n";
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
	#print "maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
			if ($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"} > 50){
				
				my %subtract_hash;
				
## used up to 4.02, when i suddenly realized that subtract_hash construction looks very suspicious				
#				if ($subtract_tallest){
#					my $tallest_tip = ${$static_hash_of_nodes{$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}}};	
#					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
#					my $index_of_ancestor_node;
#					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
#				
#					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
#						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
#							$index_of_ancestor_node = $n;
#							last;
#						}
#						my $depth = $static_distance_hash{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
#						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
#						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
#					}
#				}
##
				
				if ($subtract_tallest){
					my $tallest_tip = ${$static_hash_of_nodes{$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node; 
					my $path_length = $tallest_tip->get_branch_length; 
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last; # even if ancestors (from get_ancestors) contain the node itself, it won't make any difference
						}
						my $depth = $static_distance_hash{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				
				foreach my $bin (sort {$a <=> $b} keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					print "bin $bin adding ".$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]."\n";
					$total_muts += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($subtract_tallest && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0){
					#print "total muts $total_muts \n";
					print FILE "site $ind node ".$node->get_name()." maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
					my $check_local_lengths_sum;
					my $check_total_obs;
					my $check_total_exp;
					foreach my $bin (sort {$a <=> $b} (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}})){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						
						
						## on 4.02 changed this..
						#if ($total_length > 0){ # 30.10 - the line up there was a mistake.
						#	my $local_length = $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
						#	if ($subtract_tallest && $subtract_hash{$bin}){
						#		$local_length -= $subtract_hash{$bin};
						#	}
						#	$check_local_lengths_sum += $local_length;
						#	$hist{$site_node}{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
						#	$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
						#}
						##
						
						## ..to this
						
							my $local_length = $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($subtract_tallest && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$check_local_lengths_sum += $local_length;
							$hist{$site_node}{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
						print "adding to obs bin $bin ".$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]."\n";
						##
						
						
						##on 4.02 commented out 
						    if (!$hist{$site_node}{$bin}[0]){
						    	$hist{$site_node}{$bin}[0] += 0;
						    }
						    if (!$hist{$site_node}{$bin}[1]){
						    	$hist{$site_node}{$bin}[1] += 0;
						    }
						##
	 print FILE "$bin,".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";
	 $check_total_obs += $hist{$site_node}{$bin}[0];
	 $check_total_exp += $hist{$site_node}{$bin}[1];
				}

				if ($total_length == $check_local_lengths_sum){
				print "local lengths sumtest ok: $total_length $check_local_lengths_sum\n";
				}
				else {
				print "local length sumtest falied! total $total_length, local_sum $check_local_lengths_sum\n";
				}
				if ($check_total_obs-$check_total_exp < 0.001 && -$check_total_obs+$check_total_exp < 0.001 ){
				print "obsexp sumtest ok\n";
				}
				else {
				print "obsexp sumtest falied! total obs $check_total_obs, total exp $check_total_exp total_muts $total_muts site $ind node ".$node->get_name()."\n";
				}
				}
				else {
				print "total length 0\n";
				}
			}
			
		}
	}	
	close FILE;
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}

## Since 5.02
sub depth_groups_entrenchment_optimized_selection_alldepths {
	my $step = $_[0];
	my $ancestral_nodes = $_[2];
	my $subtract_tallest = $_[3];
	my $tag = $_[4];
	
	my $filename = $static_protein."_for_enrichment_".$tag;
	open FILE, ">>$filename";
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print FILE ">new iteration\n";
	
	my @group;
	if ($_[5]){
		@group = @{$_[5]};
	}
	else {
		@group = (1..565);
	}
	
	foreach my $key (keys %{$ancestral_nodes}){
	#print " ancestral node $key \n";
	}
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root, $subtract_tallest);
	$self->my_visit_depth_first ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);

	foreach my $ind (@group){
	#print "here my ind $ind\n";
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $nodename = $node->get_name();
			my $site_node = $ind."_".$nodename;
			#print "here my site_node $site_node\n";
			if (!$ancestral_nodes->{$nodename}){
				#print "$ind $nodename NOT IN REAL ANCESTORS\n";
				my $totmut;
				my $totlen;
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$totmut += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$totlen += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
				}	
			#	print " totmut $totmut, totlen $totlen\n";
				next;
			}
		#	print "$site_node is still here\n";
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
	#print "maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
			if ($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"} > 50){
				
				my %subtract_hash;
				
				
				if ($subtract_tallest){
					my $tallest_tip = ${$static_hash_of_nodes{$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node; 
					my $path_length = $tallest_tip->get_branch_length; 
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last; # even if ancestors (from get_ancestors) contain the node itself, it won't make any difference
						}
						my $depth = $static_distance_hash{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				
				foreach my $bin (sort {$a <=> $b} keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					#print "bin $bin adding ".$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]."\n";
					$total_muts += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($subtract_tallest && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0){
					#print "total muts $total_muts \n";
					print FILE "site $ind node ".$node->get_name()." maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
					my $check_local_lengths_sum;
					my $check_total_obs;
					my $check_total_exp;
					foreach my $bin (sort {$a <=> $b} (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}})){
						
							my $local_length = $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($subtract_tallest && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$check_local_lengths_sum += $local_length;
							$hist{$site_node}{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
	#print "adding to obs bin $bin ".$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]."\n";

						    if (!$hist{$site_node}{$bin}[0]){
						    	$hist{$site_node}{$bin}[0] += 0;
						    }
						    if (!$hist{$site_node}{$bin}[1]){
						    	$hist{$site_node}{$bin}[1] += 0;
						    }

	 print FILE "$bin,".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";
	 $check_total_obs += $hist{$site_node}{$bin}[0];
	 $check_total_exp += $hist{$site_node}{$bin}[1];
				}

				#if ($total_length == $check_local_lengths_sum){
				#print "local lengths sumtest ok: $total_length $check_local_lengths_sum\n";
				#}
				#else {
				#print "local length sumtest failed! total $total_length, local_sum $check_local_lengths_sum\n";
				#}
				#if ($check_total_obs-$check_total_exp < 0.001 && -$check_total_obs+$check_total_exp < 0.001 ){
				#print "obsexp sumtest ok\n";
				#}
				#else {
				#print "obsexp sumtest failed! total obs $check_total_obs, total exp $check_total_exp total_muts $total_muts site $ind node ".$node->get_name()."\n";
				#}
				}

			}
			
		}
	}	
	close FILE;
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}


# 3.11 prints maxdepth for each site_node - so that the necessary maxdepth could be chosen afterwards; 
# does not print site_node for each bin.

sub depth_groups_entrenchment_optimized_selector_alldepths {
	my $step = $_[0];
	my $restriction = $_[1];
	my $subtract_tallest = $_[2];

	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "real data\n";
	
	my @group;
	if ($_[3]){
		@group = @{$_[3]};
	}
	else {
		@group = (1..565);
	}
	

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root, $subtract_tallest);
	$self->my_visit_depth_first ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $site_node = $node->get_name();
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				
				my %subtract_hash;
				
				if ($subtract_tallest){
					#print "just checking ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}."\n";
					my $tallest_tip = ${$static_hash_of_nodes{$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = $static_distance_hash{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($subtract_tallest && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0 && $total_muts > 0){ # 21.10 added total_muts > 0
				#print "total muts $total_muts \n";
				print "site $ind node ".$node->get_name()." maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
					foreach my $bin (sort {$a <=> $b} (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}})){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						 if ($total_length > 0){ 
							my $local_length = $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($subtract_tallest && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$site_node}{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
						}
						if (!$hist{$site_node}{$bin}[0]){
							$hist{$site_node}{$bin}[0] += 0;
						}
						if (!$hist{$site_node}{$bin}[1]){
							$hist{$site_node}{$bin}[1] += 0;
						}
	 print "$bin,".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";
				}
				}
			}
			
		}
	}	
	
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}

# 21.10 for permutations select only ancestral mutations (with > 0 number of reversions) present in real data 

sub depth_groups_entrenchment_optimized_selection {
	my $step = $_[0];
	my $restriction = $_[1];
	my $ancestral_nodes = $_[2];
	my $subtract_tallest = $_[3];

	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	
	my @group;
	if ($_[4]){
		@group = @{$_[4]};
	}
	else {
		@group = (1..565);
	}
	
	

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root, $subtract_tallest);
	$self->my_visit_depth_first ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $site_node = $node->get_name();
			if (!$ancestral_nodes->{$site_node}){
			#	print "$ind $site_node NOT IN REAL ANCESTORS\n";
				my $totmut;
				my $totlen;
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$totmut += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$totlen += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
				}	
			#	print " totmut $totmut, totlen $totlen\n";
				next;
			}
		#	print "$site_node is still here\n";
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				
				my %subtract_hash;
				
				if ($subtract_tallest){
					my $tallest_tip = ${$static_hash_of_nodes{$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = $static_distance_hash{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($subtract_tallest && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0){
					#print "total muts $total_muts \n";
					foreach my $bin (sort {$a <=> $b} (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}})){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						if ($total_length > 0){ # 30.10 - the line up there was a mistake.
							my $local_length = $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($subtract_tallest && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$site_node}{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
						}
						if (!$hist{$site_node}{$bin}[0]){
							$hist{$site_node}{$bin}[0] += 0;
						}
						if (!$hist{$site_node}{$bin}[1]){
							$hist{$site_node}{$bin}[1] += 0;
						}
	 print "$bin,$ind,".$node->get_name().",".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";
				}
				}
			}
			
		}
	}	
	
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}

# 21.10 select meaningful ancestral mutations from real data
sub depth_groups_entrenchment_optimized_selector {
	my $step = $_[0];
	my $restriction = $_[1];
	my $subtract_tallest = $_[2];

	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	
	my @group;
	if ($_[3]){
		@group = @{$_[3]};
	}
	else {
		@group = (1..565);
	}
	

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root, $subtract_tallest);
	$self->my_visit_depth_first ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $site_node = $node->get_name();
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_subtree_info{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				print " passed 150 limit\n";
				my %subtract_hash;
				
				if ($subtract_tallest){
					#print "just checking ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}."\n";
					my $tallest_tip = ${$static_hash_of_nodes{$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = $static_distance_hash{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				foreach my $bin (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					print "added ".$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]." to muts\n";
					$total_length += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					print "added ".$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1]." to length\n";
					if ($subtract_tallest && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0 && $total_muts > 0){ # 21.10 added total_muts > 0
				print "total muts $total_muts \n";
					foreach my $bin (sort {$a <=> $b} (keys %{$static_subtree_info{$node->get_name()}{$ind}{"hash"}})){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						 if ($total_length > 0){ 
							my $local_length = $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($subtract_tallest && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$site_node}{$bin}[0] += $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected	
						}
						if (!$hist{$site_node}{$bin}[0]){
							$hist{$site_node}{$bin}[0] += 0;
						}
						if (!$hist{$site_node}{$bin}[1]){
							$hist{$site_node}{$bin}[1] += 0;
						}
	 print "$bin,$ind,".$node->get_name().",".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";
				}
				}
			}
			
		}
	}	
	
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}

     sub collect_ring_squares {
		my $self = $_[0];
		my $node = $_[1];
		my $site_index = $_[2]->[0];
		my $step = $_[2]->[1];

		my $root = $self->{static_tree} -> get_root;
		$self->my_visit_depth_first($root, (), \&has_no_mutation, \&update_ring, \($site_index, $step));
    }
    
    sub depth_groups_entrenchment_hash {
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
			if ($static_depth_hash{$ind}{$node->get_name()} > 50){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			#if ($total_length > 0 && $total_muts/$total_length < 0.005){
			if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}{$ind}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}{$ind}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					
				}
				if (!$hist{$bin}{$ind}[0]){
					$hist{$bin}{$ind}[0] += 0;
				}
				if (!$hist{$bin}{$ind}[1]){
					$hist{$bin}{$ind}[1] += 0;
				}
				
			}
			}
			}
		}
		}

	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		my $obs_sum;
		my $exp_sum;
		foreach my $s (@group){
			$obs_sum += $hist{$bin}{$s}[0];
			$exp_sum += $hist{$bin}{$s}[1];
		}
		print " up to ".$bin*$step."\t".$obs_sum."\t".$exp_sum."\n";
	}
	return %hist;
	
}


sub depth_groups_entrenchment_bootstrap {
	my $step = $_[0];
	my $root = $static_tree-> get_root;
	my @pregroup = @{$_[1]};
	my $iterations = $_[2];
	
	# both group and complement must contain only variable sites:
	my %variable_sites;
	foreach my $s(1..565){
		if ($static_nodes_with_sub{$s}){
			$variable_sites{$s} = 1;
		}
	}
	
	my %ghash;
	my @group;
	foreach my $s(@pregroup){
		if ($variable_sites{$s}){
			$ghash{$s} = 1;
			push @group, $s;
		}
	}
	my @complement;
	foreach my $s(keys %variable_sites){
		if (!$ghash{$s}){
			push @complement, $s;
		}
	}
	my @variable_sites_array = keys %variable_sites;
	my %hist = depth_groups_entrenchment_hash($step, \@variable_sites_array);
	my @obs_values = diffdiff(\%hist, \@group, \@complement);
	 
	my $counter1;
	my $counter2;

	
	for (my $i = 0; $i < $iterations; $i++){
		
		my @boot_group = shuffle @variable_sites_array;
		my @boot_complement = splice (@boot_group, scalar @group, scalar @variable_sites_array - scalar @group);
		my @boot_values = diffdiff(\%hist, \@boot_group, \@boot_complement);
		
		if ($boot_values[4] >= $obs_values[4]){
			$counter1++;
		}
		if ($boot_values[4] <= $obs_values[4]){
			$counter2++;
		}


	}
	print "epistasis enrichment ".$counter1/$iterations."\n";
	print "environment enrichment ".$counter2/$iterations."\n";

}


sub depth_groups_entrenchment_heaps {
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
			$self->my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_depth_hash{$ind}{$node->get_name()} > 50){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			if ($total_length > 0 && $total_muts/$total_length < 0.005){
			#if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				my $heap0 = 0;
				my $heap1 = 1;
				if ($static_ring_hash{$ind}{$node->get_name()}{1}[0] > $total_muts*$static_ring_hash{$ind}{$node->get_name()}{1}[1]/$total_length) {
					$heap0 = 2;
					$heap1 = 3;
				}
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					$hist{$bin}[$heap0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					$hist{$bin}[$heap1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					
				}
				if (!$hist{$bin}[$heap0]){
					$hist{$bin}[$heap0] += 0;
				}
				if (!$hist{$bin}[$heap1]){
					$hist{$bin}[$heap1] += 0;
				}
				print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
			}
			}
			}
		}
		}
	}
	
	foreach my $bin (sort {$a <=> $b} keys %hist){
		print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\t".$hist{$bin}[2]."\t".$hist{$bin}[3]."\n";
	}
	
	
}

 	sub max_depth {
 		my $self = $_[0];
 		my $node = $_[1];
		my $site_index = $_[2]->[0];
		my $starting_node = $_[2]->[2];
		my $depth = $_[3] - $starting_node->get_branch_length;
	#	print " current node is ".$node-> get_name()."\n";
	#	print " current depth is $depth\n";
		if ($node eq $starting_node){
	#		print $node-> get_name()." equal ";
			return;
		}
		if ($starting_node -> is_terminal){
	#		print $node-> get_name()." terminal ";
			return;
		}
 		
 			if ($self->{static_depth_hash}{$site_index}{$starting_node->get_name()}){
 	#			print "\n".$static_depth_hash{$site_index}{$starting_node->get_name()}.", ".$depth."; max is ";
 				$self->{static_depth_hash}{$site_index}{$starting_node->get_name()} = max($self->{static_depth_hash}{$site_index}{$starting_node->get_name()}, $depth);
 	#			print $static_depth_hash{$site_index}{$starting_node->get_name()}."\n";
 			}
 			else {
 	#			print "\n init: $depth\n";
 				$self->{static_depth_hash}{$site_index}{$starting_node->get_name()} = $depth;
 			}
 		
 	}


#stale: based on shuffling matrices
sub entrenchment_blue_violet_bootstrap{ 
	my $self = shift;
	my $iterations = $_[0];
	my $restriction = $_[1];
	$self->set_distance_matrix();
	my %matrix = $self->incidence_matrix(); #!
	$self->print_incidence_matrix(\%matrix);
	my %obs_hash = $self->depth_groups_entrenchment_optimized(1,$restriction);
	my %real_obs_hash;
	my %real_exp_hash;
	my $norm;
	
	my $maxbin = 0;
	
	foreach my $bin(keys %obs_hash){
			$maxbin = max($bin, $maxbin);
			$norm += $obs_hash{$bin}[0];
	}
	
	%static_ring_hash = ();
	%static_subtree_info = ();
	
	my @simulated_hists;
	for (my $i = 1; $i <= $iterations; $i++){
		my @mock_mutmaps = $self ->read_incidence_matrix("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash = depth_groups_entrenchment_optimized(1,$restriction);
		push @simulated_hists, \%hash;
		
		my $sum;
		foreach my $bin(1..$maxbin){
				$sum += $hash{$bin}[0];
		}
		foreach my $bin(1..$maxbin){
			$hash{$bin}[0] = $hash{$bin}[0]*$norm/$sum;
			$hash{$bin}[1] = $hash{$bin}[1]*$norm/$sum;
		}
		
		%static_ring_hash = ();
		%static_subtree_info = ();
	}

	store \@simulated_hists, "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_1_100noneed_bv_simulation_".$restriction."_stored";
	my $arref = retrieve("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_1_100noneed_bv_simulation_".$restriction."_stored");
	
	open CSV, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_1_100noneed_bv_simulation_".$restriction.".csv";
	foreach my $bin(1..$maxbin){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..$iterations-1){
		foreach my $bin(1..$maxbin){
			print CSV ($arref->[$i]->{$bin}->[0]).",".($arref->[$i]->{$bin}->[1]).",";
		}
		print CSV"\n";
	}
	close CSV;
	
	my $file = $prot."_1_100noneed_bv_boot_median_test_".$restriction;
	open FILE, ">$file";

	foreach my $bin(1..$maxbin){
			$real_obs_hash{$bin} = $obs_hash{$bin}->[0];
			$real_exp_hash{$bin} = $obs_hash{$bin}->[1];
		}
	my $real_obs_median = hist_median_for_hash(\%real_obs_hash);
	my $real_exp_median = hist_median_for_hash(\%real_exp_hash);
	
	print FILE "\n observed medians: $real_obs_median , $real_exp_median\n";
	
	my $pval_epi;
	my $pval_env;
	
	foreach my $i(0..$iterations-1){
		my %boot_obs_hash;
		my %boot_exp_hash;
		foreach my $bin(1..$maxbin){
			$boot_obs_hash{$bin} = $arref->[$i]->{$bin}->[0];
			$boot_exp_hash{$bin} = $arref->[$i]->{$bin}->[1];
		}
		my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash);
		my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash);
		print FILE "\n boot medians: $boot_obs_median , $boot_exp_median\n";
		if ($boot_obs_median - $boot_exp_median >= $real_obs_median - $real_exp_median){
			$pval_env += 1;
		}
		if ($boot_obs_median - $boot_exp_median <= $real_obs_median - $real_exp_median){
			$pval_epi += 1;
		}

	}
	print FILE "pvalue epistasis ".($pval_epi/$iterations)." pvalue environment ".($pval_env/$iterations);
	close FILE;
}

#stale: based on matrices shuffling
sub check_entrenchment_blue_violet_bootstrap {
	my $prot = $_[0];
	my $iterations = $_[1];
	set_mutmap($prot, "nsyn");
	$self->set_distance_matrix();
	my %matrix = $self->incidence_matrix(); #!
	$self -> print_incidence_matrix(\%matrix);
	my %obs_hash = depth_groups_entrenchment_optimized(1,0);
	my %real_obs_hash;
	my %real_exp_hash;
	my $norm;
	
	my $maxbin = 0;
	
	foreach my $bin(keys %obs_hash){
			$maxbin = max($bin, $maxbin);
			$norm += $obs_hash{$bin}[0];
	}
	
	%static_ring_hash = ();
	%static_subtree_info = ();
	
	my @simulated_hists;
	my %simulated_medians;
	my %all_simulated_medians;
	my %all_hash;
	
	open STORAGE, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_todelete_2434node.csv";
	for (my $i = 1; $i <= $iterations; $i++){
		my @mock_mutmaps = $self ->read_incidence_matrix("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};



		my %prehash = check_depth_groups_entrenchment_optimized(1);
		push @simulated_hists, \%prehash;
		my %hash;
		my $sum;
		
		
		foreach my $bin(1..$maxbin){
			print "----\n";
				foreach my $mutnum (keys %prehash){
					if ($mutnum > 0){
						foreach my $site_node(keys %{$prehash{$mutnum}}){
							if ($site_node eq "INTNODE2434"){
								$sum += $prehash{$mutnum}{$site_node}{$bin}[0];
								$hash{$mutnum}{"counter"} = $prehash{$mutnum}{$site_node}{"counter"};
								$hash{$mutnum}{$bin}[1] = $prehash{$mutnum}{$site_node}{$bin}[1];
								$hash{$mutnum}{$bin}[0] = $prehash{$mutnum}{$site_node}{$bin}[0];
								print "site node ".$site_node." mutnum $mutnum\n";
							}
						}
					}	
			}
		}
		if ($sum == 0){
			next;
		}
				foreach my $mutnum (keys %hash){
					if ($mutnum > 0){
					#print STORAGE ">$mutnum,".$hash{$mutnum}{"counter"}."\n";
					foreach my $bin(1..$maxbin){
						$hash{$mutnum}{$bin}[0] = $hash{$mutnum}{$bin}[0]*$norm/$sum;
						$all_hash{$mutnum}{$bin}[0] += $hash{$mutnum}{$bin}[0];
						$hash{$mutnum}{$bin}[1] = $hash{$mutnum}{$bin}[1]*$norm/$sum;
						$all_hash{$mutnum}{$bin}[1] += $hash{$mutnum}{$bin}[1];
					#	print STORAGE "$bin,".$hash{$mutnum}{$bin}[0].",".$hash{$mutnum}{$bin}[1]."\n";
					}
					$all_hash{$mutnum}{"counter"} += $hash{$mutnum}{"counter"};
					my $boot_obs_median = hist_median_for_hash_arr(\%{$hash{$mutnum}}, 0);
					push @{$simulated_medians{$mutnum}[0]}, $boot_obs_median;
					my $boot_exp_median = hist_median_for_hash_arr(\%{$hash{$mutnum}}, 1);
					push @{$simulated_medians{$mutnum}[1]}, $boot_exp_median;
					}
			}
		
		%static_ring_hash = ();
		%static_subtree_info = ();
	}
	
	my @sorted_keys = sort {$a <=> $b} keys %simulated_medians;
	my %total_medians;
	foreach my $mutnum(@sorted_keys){
		my $stat_obs = Statistics::Descriptive::Full->new();
		$stat_obs->add_data(\@{$simulated_medians{$mutnum}[0]});
		my $stat_exp = Statistics::Descriptive::Full->new();
		$stat_exp->add_data(\@{$simulated_medians{$mutnum}[1]});
		print STORAGE $mutnum.",".hist_median_for_hash_arr(\%{$all_hash{$mutnum}}, 0).",".
		hist_mean_for_hash_arr(\%{$all_hash{$mutnum}}, 0).",".stddev(\@{$simulated_medians{$mutnum}[0]}).",".
		$stat_obs->percentile(5).",".
		$stat_obs->percentile(95).",".
		hist_median_for_hash_arr(\%{$all_hash{$mutnum}}, 1).",".
		hist_mean_for_hash_arr(\%{$all_hash{$mutnum}}, 1).",".stddev(\@{$simulated_medians{$mutnum}[1]}).",".
		$stat_exp->percentile(5).",".
		$stat_exp->percentile(95).",".
		$all_hash{$mutnum}{"counter"}."\n";
	}
	
	
	
#	foreach my $mutnum(@sorted_keys){
	#	print STORAGE $mutnum."\n";
	#	for (my $i = 1; $i <= $iterations; $i++){
	#		print STORAGE $simulated_medians{$mutnum}[0][$i].",".$simulated_medians{$mutnum}[1][$i]."\n";
	#	}
	#}
close STORAGE;
	#store \@simulated_hists, "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_check_simulation_150_hash_stored";
	
}

# active on 26.10
#check_entrenchment_blue_violet_bootstrap ("h1", 50);

#	check_part_2("h3", 20, 150);


            
            


sub check_part_2	{
		my $prot = $_[0];
	my $iterations = $_[1];
	my $restriction = $_[2];

	my $maxbin = 400;
		my $arref = retrieve("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_check_simulation_".$restriction."_stored");
	
	my $file = $prot."_check_statistics_".$restriction;
	open FILE, ">$file";
	
	my $pval_epi;
	my $pval_env;
	foreach my $mutnum (1..50){
		print $mutnum."\n";
			my %complete_boot_obs_hash;
			my %complete_boot_exp_hash;
			my $complete_counter;
		foreach my $i(0..$iterations-1){
			my %boot_obs_hash;
			my %boot_exp_hash;
			
			my $counter;
			foreach my $maxdepth (0..7){
				#print FILE "\t maxdepth $maxdepth\n";
				if ($maxdepth >= $restriction/50){
					#print FILE "maxdepth ok\n";
					if ($arref->[$i]->[$maxdepth]){
						foreach my $testkey (keys %{$arref->[$i]->[$maxdepth]}){
							#print " found mutnum $testkey\n";
						}
						foreach my $site_node (keys %{$arref->[$i]->[$maxdepth]->{$mutnum}}){
							#print FILE " \t site_node $site_node\n";
								$counter++;
								$complete_counter++;
								foreach my $bin(1..$maxbin){
									$boot_obs_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[0];
									$boot_exp_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[1];
									$complete_boot_obs_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[0];
									$complete_boot_exp_hash{$bin} += $arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[1];
									#print FILE " added".$arref->[$i]->[$maxdepth]->{$mutnum}->{$site_node}->{$bin}->[0]." and ".$arref->[$i]->{$maxdepth}->{$mutnum}->{$site_node}->{$bin}->[1]."\n";
								
								}
						
					}
					}
				}
			}
			
			
			my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash);
			my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash);
			
			#print FILE "\n boot medians: $boot_obs_median , $boot_exp_median".."\n";

		}
		my $complete_boot_obs_median = hist_median_for_hash(\%complete_boot_obs_hash);
		my $complete_boot_exp_median = hist_median_for_hash(\%complete_boot_exp_hash);
		print FILE "mutnum $mutnum obs median $complete_boot_obs_median, exp median $complete_boot_exp_median, counter $complete_counter\n";
	}

	close FILE;
}


sub entrenchment_bootstrap_full_selection_vector { 
	my $self = shift;
	my $iterations = shift;
	my $restriction = shift;
	$self->set_distance_matrix(); # reads file from data folder
	my %matrix = $self->incidence_matrix(); #!
	$self -> print_incidence_matrix(\%matrix);
	my %obs_hash = $self->depth_groups_entrenchment_optimized_selector_alldepths_2(1,$restriction); #bin, restriction #09.08 changed depth_groups_entrenchment_optimized_selector_alldepths (which was deleted) to depth_groups_entrenchment_optimized_selector_alldepths_2
	my %ancestor_nodes;
	foreach my $ancnode(keys %obs_hash){
		$ancestor_nodes{$ancnode} = 1;
	#	print "ANCNODE ".$ancnode."\n";
	}
	
	my $norm;
	my $maxbin = 0;
	
	foreach my $site_node(keys %obs_hash){
		foreach my $bin(keys %{$obs_hash{$site_node}}){
			$maxbin = max($bin, $maxbin);
			$norm += $obs_hash{$site_node}{$bin}[0];
		}
	}

	my $mock_mutmap = myclone();  
	my @simulated_hists;
	
	for (my $i = 1; $i <= $iterations; $i++){
		$mock_mutmap->shuffle_mutator();
		my %hash;
		my $sum;
		my %prehash = $mock_mutmap->depth_groups_entrenchment_optimized_selection_alldepths(1,$restriction,\%ancestor_nodes, "overwrite"); #bin, restriction, ancestor_nodes
		
		foreach my $bin(1..$maxbin){
				foreach my $site_node(keys %prehash){
					#print " READING $site_node\n";
								$sum += $prehash{$site_node}{$bin}[0];
								$hash{$bin}[1] += $prehash{$site_node}{$bin}[1];
								$hash{$bin}[0] += $prehash{$site_node}{$bin}[0];
								#print "site node ".$site_node." mutnum $mutnum\n";			
			}
		}
		foreach my $bin(1..$maxbin){
			$hash{$bin}[0] = $hash{$bin}[0]*$norm/$sum;
			$hash{$bin}[1] = $hash{$bin}[1]*$norm/$sum;
		}
		push @simulated_hists, \%hash;
	}
	
	store \@simulated_hists, $self->{static_output_base}.$self->{static_protein}."_testselector_vector_".$restriction."_stored";
	my $arref = retrieve($self->{static_output_base}.$self->{static_protein}."_testselector_vector_".$restriction."_stored");
	my $csvfile = $self->{static_output_base}.$self->{static_protein}."_testselector_vector_".$restriction.".csv";
	open CSV, ">$csvfile";
	foreach my $bin(1..$maxbin){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..$iterations-1){
		foreach my $bin(1..$maxbin){
			print CSV ($arref->[$i]->{$bin}->[0]).",".($arref->[$i]->{$bin}->[1]).",";
		}
		print CSV"\n";
	}
	close CSV;
	
	my $testfile = $self->{static_output_base}.$self->{static_protein}."_testselector_vector_boot_median_test_".$restriction;
	open FILE, ">$testfile";
	my %flat_obs_hash;
	foreach my $bin(1..$maxbin){
		foreach my $node (keys %obs_hash){
			$flat_obs_hash{$bin} += $obs_hash{$node}{$bin}[0];
		}
	}
	my $obs_median = hist_median_for_hash(\%flat_obs_hash);
	
	print FILE "\n observed median: $obs_median\n";
	
	my $pval_epi;
	my $pval_env;
	
	foreach my $i(0..$iterations-1){
		my %hash;
		foreach my $bin(1..$maxbin){
			$hash{$bin} = $arref->[$i]->{$bin}->[0];
		}
		my $boot_median = hist_median_for_hash(\%hash);
		print FILE "\n boot median: $boot_median\n";
		if ($boot_median <= $obs_median){
			$pval_epi += 1;
		}
		if ($boot_median >= $obs_median){
			$pval_env += 1;
		}

	}
	print FILE "pvalue epistasis ".($pval_epi/$iterations)." pvalue environment ".($pval_env/$iterations);
	close FILE;
}


sub enrichment_blue_violet_optimized {
	my $obsfile = $_[0];
	my $bootfile = $_[1];
	my @group = @{$_[2]};
	
	my %grouphash;
	foreach my $ind(@group){
		$grouphash{$ind} = 1;
	}
	
	my %boot_obshash;
	my %boot_exphash;
	my %obshash;
	my %exphash;

	
	open FILE, "<$obsfile" or die "Cannot open $obsfile\n";
	
	my $maxbin = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		
		$maxbin = max($bin, $maxbin);
		# bin, ind, node, obs, exp
		$obshash{$ind}{$bin} += $arr[3];
		$exphash{$ind}{$bin} += $arr[4];
	}

	my $norm;
	foreach my $bin(1..$maxbin){
		foreach my $ind (keys %obshash){
			$norm += $obshash{$ind}{$bin};
		}
	}
	
	my %obshash_g;
	my %obshash_c;
	my %exphash_g;
	my %exphash_c;
	
	my $counter_g;
	my $counter_c;
	
	foreach my $bin(1..$maxbin){
			foreach my $ind(keys %obshash){
				if ($grouphash{$ind}){
					$obshash_g{$bin} += $obshash{$ind}{$bin}; # group
					$exphash_g{$bin} += $exphash{$ind}{$bin}; # group
				}
				else {
					$obshash_c{$bin} += $obshash{$ind}{$bin}; # complement
					$exphash_c{$bin} += $exphash{$ind}{$bin}; # complement
				}
			}
	}
		
		
		foreach my $ind(keys %obshash){
				if ($grouphash{$ind}){
					if ($obshash{$ind}){
						$counter_g ++;
					}
				}
				else {
					if ($obshash{$ind}){
						$counter_c ++;
					}
				}
			}
		
		
	print "number of sites: group $counter_g , complement $counter_c\n";
		
	my $obs_median_g = hist_median_for_hash(\%obshash_g);
	my $obs_median_c = hist_median_for_hash(\%obshash_c);
	my $exp_median_g = hist_median_for_hash(\%exphash_g);
	my $exp_median_c = hist_median_for_hash(\%exphash_c);
	print " obs_median_g $obs_median_g  obs_median_c $obs_median_c exp_median_g $exp_median_g  exp_median_c $exp_median_c\n ";	
	
	open FILE, "<$bootfile" or die "Cannot open $bootfile\n";
	my $iteration = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			$iteration++;
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		# bin, ind, node, obs,exp
		$boot_obshash{$iteration}{$ind}{$bin} += $arr[3];
		$boot_exphash{$iteration}{$ind}{$bin} += $arr[4];
	}
	close FILE;
	print "$iteration iterations\n";
	my $pval_env_enrichment;
	my $pval_env_depletion;
	my $pval_epi_enrichment;
	my $pval_epi_depletion;
	

	
	
	foreach my $i(keys %boot_obshash){
		my $sum;
		my %boot_obshash_g;
		my %boot_obshash_c;
		my %boot_exphash_g;
		my %boot_exphash_c;
		
		my $counter_g;
		my $counter_c;
			foreach my $ind(keys %{ $boot_obshash{$i}}){
				if ($grouphash{$ind}){
					if ($boot_obshash{$i}{$ind}){
						$counter_g ++;
					}
				}
				else {
					if ($boot_obshash{$i}{$ind}){
						$counter_c ++;
					}
				}
			}
		
		
	print "number of sites: group $counter_g , complement $counter_c\n";
		
		
		
		foreach my $bin(1..$maxbin){
			foreach my $ind (1..565){
				$sum += $boot_obshash{$i}{$ind}{$bin};
			}
		}
		foreach my $bin(1..$maxbin){
			foreach my $ind(keys %{$boot_obshash{$i}}){
				if ($grouphash{$ind}){
					$boot_obshash_g{$bin} += $boot_obshash{$i}{$ind}{$bin}*$norm/$sum; # group
					$boot_exphash_g{$bin} += $boot_exphash{$i}{$ind}{$bin}*$norm/$sum; # group
				}
				else {
					$boot_obshash_c{$bin} += $boot_obshash{$i}{$ind}{$bin}*$norm/$sum; # complement
					$boot_exphash_c{$bin} += $boot_exphash{$i}{$ind}{$bin}*$norm/$sum; # complement
				}
			}
		}
		my $boot_obs_median_g = hist_median_for_hash(\%boot_obshash_g);
		my $boot_obs_median_c = hist_median_for_hash(\%boot_obshash_c);
		my $boot_exp_median_g = hist_median_for_hash(\%boot_exphash_g);
		my $boot_exp_median_c = hist_median_for_hash(\%boot_exphash_c);
		print " boot_obs_median_g $boot_obs_median_g  boot_obs_median_c $boot_obs_median_c boot_exp_median_g $boot_exp_median_g  boot_exp_median_c $boot_exp_median_c\n ";	
		
		if (($boot_obs_median_g-$boot_exp_median_g) - ($boot_obs_median_c-$boot_exp_median_c) 
			 >= ($obs_median_g-$exp_median_g) - ($obs_median_c-$exp_median_c) ){
			$pval_env_enrichment += 1;
		}
		if (($boot_obs_median_g-$boot_exp_median_g) - ($boot_obs_median_c-$boot_exp_median_c) 
			 <= ($obs_median_g-$exp_median_g) - ($obs_median_c-$exp_median_c)){
			$pval_env_depletion += 1;
		}
		if (-($boot_obs_median_g-$boot_exp_median_g) + ($boot_obs_median_c-$boot_exp_median_c) 
			 >= -($obs_median_g-$exp_median_g) + ($obs_median_c-$exp_median_c)){
			$pval_epi_enrichment += 1;
		}
		if (-($boot_obs_median_g-$boot_exp_median_g) + ($boot_obs_median_c-$boot_exp_median_c) 
			 <= -($obs_median_g-$exp_median_g) + ($obs_median_c-$exp_median_c)){
			$pval_epi_depletion += 1;
		}
		
	}
	
	print " pval epi enrichment ".$pval_epi_enrichment/$iteration."\n";
	print " pval epi depletion ".$pval_epi_depletion/$iteration."\n";
	print " pval env enrichment ".$pval_env_enrichment/$iteration."\n";
	print " pval env depletion ".$pval_env_depletion/$iteration."\n";
}


sub enrichment_optimized {
	my $obsfile = $_[0];
	my $bootfile = $_[1];
	my @group = @{$_[2]};
	
	my %grouphash;
	foreach my $ind(@group){
		$grouphash{$ind} = 1;
	}
	
	my %hash;
	my %obshash;

	
	open FILE, "<$obsfile" or die "Cannot open $obsfile\n";
	my $maxbin = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		$maxbin = max($bin, $maxbin);
		# bin, ind, node, obs
		$obshash{$ind}{$bin} += $arr[3];
	}

	my $norm;
	foreach my $bin(1..$maxbin){
		foreach my $ind (keys %obshash){
			$norm += $obshash{$ind}{$bin};
		}
	}
	
	my %obshash_g;
	my %obshash_c;
	foreach my $bin(1..$maxbin){
			foreach my $ind(keys %obshash){
				if ($grouphash{$ind}){
					$obshash_g{$bin} += $obshash{$ind}{$bin}; # group
				}
				else {
					$obshash_c{$bin} += $obshash{$ind}{$bin}; # complement
				}
			}
	}
		
	my $obs_median_g = hist_median_for_hash(\%obshash_g);
	my $obs_median_c = hist_median_for_hash(\%obshash_c);
	print " obs_median_g $obs_median_g  obs_median_c $obs_median_c\n ";	
	
	open FILE, "<$bootfile" or die "Cannot open $bootfile\n";
	my $iteration = 0;
	while(<FILE>){
		if ($_ =~ /^rad/){
			$iteration++;
			next;
		}
		my @arr = split(/,/);
		my $ind = $arr[1] ;
		my $bin = $arr[0];
		# bin, ind, node, obs
		$hash{$iteration}{$ind}{$bin} += $arr[3];
	}
	close FILE;
	print "$iteration iterations\n";
	my $pval_env_enrichment;
	my $pval_env_depletion;
	my $pval_epi_enrichment;
	my $pval_epi_depletion;
	
	
	foreach my $i(keys %hash){
		my $sum;
		my %hash_g;
		my %hash_c;
		foreach my $bin(1..$maxbin){
			foreach my $ind (1..565){
				$sum += $hash{$i}{$ind}{$bin};
			}
		}
		foreach my $bin(1..$maxbin){
			foreach my $ind(keys %{$hash{$i}}){
				if ($grouphash{$ind}){
					$hash_g{$bin} += $hash{$i}{$ind}{$bin}*$norm/$sum; # group
				}
				else {
					$hash_c{$bin} += $hash{$i}{$ind}{$bin}*$norm/$sum; # complement
				}
			}
		}
		my $boot_median_g = hist_median_for_hash(\%hash_g);
		my $boot_median_c = hist_median_for_hash(\%hash_c);
		print " boot_median_g $boot_median_g  boot_median_c $boot_median_c\n ";	
		
		if ($obs_median_g-$obs_median_c <= $boot_median_g-$boot_median_c){
			$pval_env_enrichment += 1;
		}
		if ($obs_median_g-$obs_median_c >= $boot_median_g-$boot_median_c){
			$pval_env_depletion += 1;
		}
		if (-$obs_median_g+$obs_median_c <= -$boot_median_g+$boot_median_c){
			$pval_epi_enrichment += 1;
		}
		if (-$obs_median_g+$obs_median_c >= -$boot_median_g+$boot_median_c){
			$pval_epi_depletion += 1;
		}
		
	}
	
	print " pval epi enrichment ".$pval_epi_enrichment/$iteration."\n";
	print " pval epi depletion ".$pval_epi_depletion/$iteration."\n";
	print " pval env enrichment ".$pval_env_enrichment/$iteration."\n";
	print " pval env depletion ".$pval_env_depletion/$iteration."\n";
		
}

sub depth_groups_entrenchment_optimized {
	my $self = shift;
	my $step = $_[0];
	my $restriction = $_[1];
	my $root = $self->{static_tree}-> get_root;
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

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root);
	$self->visitor_coat($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
		foreach my $node(@{$self->{static_nodes_with_sub}{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				foreach my $bin (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0){
					foreach my $bin (sort {$a <=> $b} (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}})){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						if ($total_length > 0){ # 30.10 the ine up there was a mistake	
							$hist{$bin}[0] += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$bin}[1] += $total_muts*$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1]/$total_length; #expected	
						}
						if (!$hist{$bin}[0]){
							$hist{$bin}[0] += 0;
						}
						if (!$hist{$bin}[1]){
							$hist{$bin}[1] += 0;
						}
	 print "$bin,$ind,".$node->get_name().",".$hist{$bin}[0].",".$hist{$bin}[1]."\n";
				}
				}
			}
			
		}
	}	
	
	#foreach my $bin (sort {$a <=> $b} keys %hist){
	#	print " up to ".$bin*$step."\t".$hist{$bin}[0]."\t".$hist{$bin}[1]."\n";
	#}
	return %hist;	
}
