sub entrenchment_bootstrap{ 
	my $prot = $_[0];
	my $restriction = $_[1];
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %matrix = incidence_matrix(); #!
	$self -> print_incidence_matrix(\%matrix);
	my %hash = depth_groups_entrenchment_optimized(10,$restriction);
	my $norm;
	foreach my $bin(1..32){
			$norm += $hash{$bin}[0];
	}

	my @simulated_hists;
	for (my $i = 1; $i < 101; $i++){
		my @mock_mutmaps = $self ->read_incidence_matrix("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash = depth_groups_entrenchment_optimized(10,$restriction);
		push @simulated_hists, \%hash;
		%static_ring_hash = ();
		%static_depth_hash = ();
	}

	store \@simulated_hists, "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_full_simulation_".$restriction."_stored";
	my $arref = retrieve("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_full_simulation_".$restriction."_stored");
	
	open CSV, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_full_simulation_".$restriction.".csv";
	foreach my $bin(1..32){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..99){
		my $sum;
		foreach my $bin(1..32){
				$sum += $arref->[$i]->{$bin}->[0];
		}
		foreach my $bin(1..32){
			print CSV ($arref->[$i]->{$bin}->[0])*$norm/$sum.",".($arref->[$i]->{$bin}->[1])*$norm/$sum.",";
		}
		print CSV"\n";
	}
	close CSV;
}

sub entrenchment_bootstrap_full{ 
	my $prot = $_[0];
	my $iterations = $_[1];
	my $restriction = $_[2];
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %matrix = incidence_matrix(); #!
	$self -> print_incidence_matrix(\%matrix);
	my %obs_hash = depth_groups_entrenchment_optimized(10,$restriction);
	my $norm;
	foreach my $bin(1..32){
			$norm += $obs_hash{$bin}[0];
	}
	
	%static_ring_hash = ();
	%static_depth_hash = ();
	%static_subtree_info = ();
	
	my @simulated_hists;
	for (my $i = 1; $i <= $iterations; $i++){
		my @mock_mutmaps = $self ->read_incidence_matrix("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash = depth_groups_entrenchment_optimized(10,$restriction);
		push @simulated_hists, \%hash;
		
		my $sum;
		foreach my $bin(1..32){
				$sum += $hash{$bin}[0];
		}
		foreach my $bin(1..32){
			$hash{$bin}[0] = $hash{$bin}[0]*$norm/$sum;
			$hash{$bin}[1] = $hash{$bin}[1]*$norm/$sum;
		}
		
		%static_ring_hash = ();
		%static_depth_hash = ();
		%static_subtree_info = ();
	}

	store \@simulated_hists, "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_noneed_simulation_".$restriction."_stored";
	my $arref = retrieve("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_noneed_simulation_".$restriction."_stored");
	
	open CSV, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_noneed_simulation_".$restriction.".csv";
	foreach my $bin(1..32){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	foreach my $i(0..$iterations-1){
		foreach my $bin(1..32){
			print CSV ($arref->[$i]->{$bin}->[0]).",".($arref->[$i]->{$bin}->[1]).",";
		}
		print CSV"\n";
	}
	close CSV;
	
	my $file = $prot."_noneed_boot_median_test_".$restriction;
	open FILE, ">$file";
	foreach my $bin(1..32){
			$obs_hash{$bin} = $obs_hash{$bin}->[0];
		}
	my $obs_median = hist_median_for_hash(\%obs_hash);
	
	print FILE "\n observed median: $obs_median\n";
	
	my $pval_epi;
	my $pval_env;
	
	foreach my $i(0..$iterations-1){
		my %hash;
		foreach my $bin(1..32){
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

# 21.10 only ancestor mutations present in real data are selected from permutations for expectance construction

sub entrenchment_bootstrap_full_selection { 
	my $prot = $_[0];
	my $iterations = $_[1];
	my $restriction = $_[2];
	my $subtract_tallest = $_[3];
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %matrix = incidence_matrix(); #!
	$self -> print_incidence_matrix(\%matrix);
	my %obs_hash = depth_groups_entrenchment_optimized_selector_alldepths(1,$restriction,$subtract_tallest); #bin, restriction, subtract-tallest
	my %ancestor_nodes;
	foreach my $ancnode(keys %obs_hash){
		$ancestor_nodes{$ancnode} = 1;
		#print "ANCNODE ".$ancnode."\n";
	}
	
	
	my $norm;
		
	my $maxbin = 0;
	
	foreach my $site_node(keys %obs_hash){
		foreach my $bin(keys %{$obs_hash{$site_node}}){
			$maxbin = max($bin, $maxbin);
			$norm += $obs_hash{$site_node}{$bin}[0];
		}
	}

	
	%static_ring_hash = ();
	%static_depth_hash = ();
	%static_subtree_info = ();
	
	my @simulated_hists;
	for (my $i = 1; $i <= $iterations; $i++){
		my @mock_mutmaps = $self ->read_incidence_matrix("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_ShuffledMatrices/$i");
		
		%static_subs_on_node = %{$mock_mutmaps[0]};
		%static_nodes_with_sub = %{$mock_mutmaps[1]};

		my %hash;
		my $sum;
		my %prehash = depth_groups_entrenchment_optimized_selection_alldepths(1,$restriction,\%ancestor_nodes,$subtract_tallest, "tag"); #bin, restriction, ancestor_nodes, subtract-tallest
		
		
		foreach my $bin(1..$maxbin){
				foreach my $site_node(keys %prehash){
					#print " READING $site_node\n";
								$sum += $prehash{$site_node}{$bin}[0];
								#$hash{$mutnum}{"counter"} = $prehash{$mutnum}{$site_node}{"counter"};
								$hash{$bin}[1] += $prehash{$site_node}{$bin}[1];
								$hash{$bin}[0] += $prehash{$site_node}{$bin}[0];
								#print "site node ".$site_node." mutnum $mutnum\n";
							
			}
		}
		
		
		#foreach my $bin(1..$maxbin){
		#		$sum += $hash{$bin}[0];
		#}
		
		foreach my $bin(1..$maxbin){
			$hash{$bin}[0] = $hash{$bin}[0]*$norm/$sum;
			$hash{$bin}[1] = $hash{$bin}[1]*$norm/$sum;
		}
		push @simulated_hists, \%hash;
		
		%static_ring_hash = ();
		%static_depth_hash = ();
		%static_subtree_info = ();
	}

	store \@simulated_hists, "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_selector_".$restriction."_stored";
	my $arref = retrieve("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_selector_".$restriction."_stored");
	
	open CSV, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_selector_".$restriction.".csv";
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
	
	my $file = $prot."_selector_boot_median_test_".$restriction;
	open FILE, ">$file";
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