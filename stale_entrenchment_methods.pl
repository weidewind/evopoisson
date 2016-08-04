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


