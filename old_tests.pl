sub test_mutmaps{
	
	foreach my $ind(keys %static_nodes_with_sub){
		my $nodes_count = scalar @{$static_nodes_with_sub{$ind}};
		print " index $ind, $nodes_count nodes ";
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			print $$node->get_name()."\t";
		}
		print "\n";
	}
	
	foreach my $nodname(keys %static_subs_on_node){
		my $sites_count = scalar keys $static_subs_on_node{$nodname};
		print " node $nodname, $sites_count sites\n";
	}
}

sub test {
my %mm = mutmap_from_files("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1.l.r.newick","/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1.all.fa");  
print Dumper ($mm{"INTNODE2340"}[0]);
print ($mm{"INTNODE2340"}[0]->{"Substitution::position"});
}

sub test2 {
my %mm = mutmap_from_files("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1.l.r.newick","/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1.all.fa");  
print Dumper ($mm{"INTNODE2340"});
#print ($mm{"INTNODE2340"}[0]->{"Substitution::position"});
}


sub testv {
my $vec1 = Bit::Vector->new(10);
$vec1->Bit_On(3);
my $string = $vec1->to_Bin();
print "'$string'\n";
my $vect = Bit::Vector->new(11);
$vect -> Interval_Copy($vec1,0,3,10);
$vect -> Bit_On(10);
$string = $vect->to_Bin();
print "'$string'\n";
$vect->Move_Left(1);
$string = $vect->to_Bin();
print "'$string'\n";
my $vect2 = Bit::Vector->new(6);
$vect2 -> Bit_On(3);
my $vect3 = Bit::Vector->new(6);
$vect3 -> Bit_On(2);
my $vect4 = Bit::Vector->new(6);
$vect4->Divide($vect2, $vect3, Bit::Vector->new(6));
print ($vect2->to_Bin()."\t".$vect3->to_Bin()."\t".$vect4->to_Bin()."\n");

}



sub patrtest{
		my $tree = parse_tree("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/n1.l.r.newick");
	foreach my $node(@{$tree->get_internals()}){
		if ($node->get_name() eq "INTNODE2473"){
		foreach my $n(@{$tree->get_internals()}){
			if ($n->get_name() eq "INTNODE2479"){
				print get_mrcn($node, $n)->get_name()."\n";
			print $node->get_name()."\t".$n->get_name()."\t".calc_true_patristic_distance($node, $n)."\n";
			}
	}
		}
	}
}

sub testero {
	my @a1 = (2,3,1,4,5,9,9,9);
	my $m = my_median(\@a1);
	print "mymedian $m ";
	my @a2 = (1,1,1,3,5,10,5);
	my $t = median_difference(\@a1, \@a2);
	print "testero $t";
}

sub medtest {
	my @array = (5,5,5,5);
	my %hash;
	$hash{1} = 2;
	$hash{2} = 4;
	$hash{3} = 4;
	$hash{4} = 10;
	print hist_median(\@array)."\n";
	print hist_median_for_hash(\%hash);
}

#test_hist_median_group();

sub test_hist_median_group {
	my @hist = (
				[1,2,0,7,8,9],
				[0,2,1,7,9,8],
				[8,7,9,0,2,1],
	);
	my @group1= (3,4,5);
	my @t1= (24,24,3);
	my @group2= (1,0,2);
	my @t2= (3,3,24);
	my  @group3= (1,0,2, 3,5,4);
	my @t3= (27,27,27);
	print hist_median_group(\@hist, \@group1)."\t" ;
	print hist_median(\@t1)."\t";
	print hist_median_group(\@hist, \@group2)."\t" ;
	print hist_median(\@t2)."\t";
	print hist_median_group(\@hist, \@group3)."\t" ;
	print hist_median(\@t3)."\t";
}

sub test_group_bootstrap{
		my @meaningful_sites = (1,2,3,6,7,8,11,12,13);
		my @group = (1,2,3);
		my @bootstrap_group = shuffle @meaningful_sites;
		my @bootstrap_complement = splice (@bootstrap_group, scalar @group, scalar @meaningful_sites - scalar @group);
		
		foreach my $g(@bootstrap_group){
			print $g.", ";
		}
		print "\n";
		foreach my $c(@bootstrap_complement){
			print $c.", ";
		}
}


sub subtreetest{
	my $tree = parse_tree("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1.l.r.newick");
my %fasta = parse_fasta("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1.all.fa");
my @mutmaps = mutmap($tree, \%fasta);
my %subs_on_node = %{$mutmaps[0]};
my %nodes_with_sub = %{$mutmaps[1]};

 compute_all_distances_in_subtree($tree, \%subs_on_node);
 my $root = $tree->get_root();
	$root->visit_depth_first(
				-post => sub{ #all daughters have been processed
					my $node=shift;
					my %dists = compute_all_distances_global($node, \%subs_on_node);
					print "node name ";
					print $node->get_name;
					print "\n";
					for my $site(keys %dists){
						print ("site name ".$site);
						print "\n";
						for my $d(@{$dists{$site}}){
							print ($d."_");
							print "\t";
						}
					}
					print "\n";
				}
	);
}

# 	test_max_depth();
 	sub test_max_depth{
 		foreach my $ind (3..3){
 		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			$node = ${$node};
			if (!$node->is_terminal){
				my @array;
				my @args = ($ind, 1, $node);
				print "Now measuring max depth for ".$node->get_name()."\n";
				my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
				print $node->get_name()."\t".$static_depth_hash{$ind}{$node->get_name()}."\n";
			}
 		}
 		}
 	}
 	
 	#test_my_visit_depth_first();

sub test_tree_lengths {
	my $tree = parse_tree("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/Mock/h1.l.r.newick");
	my $protein_name = "h1";
	tree_lengths($protein_name,$tree);
}

sub test_subtree_lengths {
	my $tree = parse_tree("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/Mock/h1.l.r.newick");
	my $protein_name = "h1";
	my $hash = subtree_lengths($protein_name, $tree);
	print $hash->{"alastar7"}." must be 13";
}

sub test_my_visit_depth_first {
	    set_mutmap("h1", "nsyn"); # for this sub to work properly path in set_mutmap had to be changed to Mock
		my $tree = parse_tree("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/Mock/h1.l.r.newick");
	    my @array;
	    my $root = $tree-> get_root;
	    my @args = (2,5, $root);
	    $self->my_visit_depth_first ($root, \@array,\&update_ring,\&has_no_mutation,\@args,0);
	    print "1: up to 5 ".$static_ring_hash{2}{"alister1"}{1}[0]."\t_".$static_ring_hash{2}{"alister1"}{1}[1]."\n";
	    print "2: up to 10 ".$static_ring_hash{2}{"alister1"}{2}[0]."\t_".$static_ring_hash{2}{"alister1"}{2}[1]."\n";
	    
}

#bullshit_test();
sub bullshit_test{
	
	set_mutmap("h1", "nsyn");
	my %matrix = incidence_matrix();
	$self -> print_incidence_matrix(\%matrix);
	my $step = 2;
	my $root = $static_tree-> get_root;
	my @array;
	my %hist;
	print "radius,site,node,observed,expected\n";
	my @group;

	foreach my $ind (1..7){
		print "\nindex $ind\n";
		foreach my $node(@{$static_nodes_with_sub{$ind}}){
			
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			print "\n   node ".$node->get_name()."\n";
			if (!$node->is_terminal){
			my @args = ($ind, $step, $node);
			$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);
			$self->my_visit_depth_first ($node, \@array,\&max_depth,\&has_no_mutation,\@args,0);
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($static_depth_hash{$ind}{$node->get_name()} > 0){
			foreach my $bin (keys %{$static_ring_hash{$ind}{$node->get_name()}}){
				print "         bin ".$bin*$step."\n";
				print "            adding ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." to muts\n";
				print "            adding ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]." to length\n";
				$total_muts += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0];
				$total_length += $static_ring_hash{$ind}{$node->get_name()}{$bin}[1];
			}
			#if ($total_length > 0 && $total_muts/$total_length < 0.005){
			if ($total_length > 0){
			foreach my $bin (sort {$a <=> $b} (keys %{$static_ring_hash{$ind}{$node->get_name()}})){
				
				if ($total_length > 0 && $static_ring_hash{$ind}{$node->get_name()}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
					print "         bin ".$bin*$step."\n";
					$hist{$bin}[0] += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]; #observed
					print "            adding ".$static_ring_hash{$ind}{$node->get_name()}{$bin}[0]." to observed\n";
					$hist{$bin}[1] += $total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length; #expected
					print "            adding ".$total_muts*$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]/$total_length." to expected\n";
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
}

## For checking iterations_gulp

sub test_gulp {
	my $fullpath = $_[0];
	my $md = $_[1];
	open GULP, "<$fullpath" or die "Cannot open $fullpath";
	
	my $site; # careful	
	my $node_name;
		while (<GULP>){

			my $max_depth;
			if ($_ =~ /^>/){
				#$iteration_number++;
		#		my $str = <GULP>;
				#print "str ".$str."\n";
		#		my @str_array = split(/\s+/, $str);
		#		my $site = $str_array[1];
		#		my $node_name = $str_array[3];
		#		$max_depth = $str_array[5];
				
				#print $site." site,".$node_name." node,".$max_depth." maxdepth\n";
			}
			my @str_array;
			my $str = <GULP>; 

			
			my %sums;
			my %hash;
			
			my $test_obs_summ;
			my $test_exp_summ;
			
			while ($str =~ /^[^>]/){ 

			
				if ($str =~ /^site/){
				
			if ($test_obs_summ == $test_exp_summ){
				print "summtest ok\n";
			}
			else {
				print "summtest failed! $site $node_name obssum $test_obs_summ, expsum $test_exp_summ\n";
			}
			$test_obs_summ = 0;
			$test_exp_summ = 0;	
				
					my @str_array = split(/\s+/, $str);
					$site = $str_array[1]; # careful
					$node_name = $str_array[3];
					$max_depth = $str_array[5];
					#print $site." site,".$node_name." node,".$max_depth." maxdepth\n";
					#print "test str $str\n";
					$str = <GULP>; ## What?! change this line to my $str = <GULP>, and one line will be skipped
					#print "test2 str $str\n";
				}
				@str_array = split(/,/, $str);
				
								#print " Maxdepth $max_depth itnum $iteration_number bin ".$str_array[0]." exp ".$str_array[2]." obs ".$str_array[1]." \n";
								
									if ($max_depth > $md){
										#foreach my $group_number(0..scalar @groups-1){
											#if ($group_hashes{$md}[$group_number]{$node_name}){
											#print "group number $group_number md $md node name $node_name\n";
												#$sums{$md}[$group_number] += $str_array[1];
												#$hash{$md}[$group_number]{$str_array[0]}[1] += $str_array[2];
												#$hash{$md}[$group_number]{$str_array[0]}[0] += $str_array[1];
												$test_obs_summ += $str_array[1];
												$test_exp_summ += $str_array[2];
												#print $hash{$md}[$group_number]{$str_array[0]}[0]." obs integral\n";
											#}
										#}
									}
								

				
				$str = <GULP>;
				#print "test3 str $str\n";
				
			}
			

			# maxbins are different for every iteration. Find maximum and use it.

			#$maxbin = max($maxbin, $str_array[0]);
#print "maxbin $maxbin, iteration number $iteration_number\n";	
#print "sum50 $sum50 sum100 $sum100 sum150 $sum150 norm 50 $norm50 norm 100 $norm100 norm 150 $norm150\n";		
			
			#foreach my $md(@maxdepths){ 
				#foreach my $group_number(0..scalar @groups-1){
					#print "maxdepth $md group number $group_number \n";
					#if ($sums{$md}[$group_number] == 0){
						#foreach my $bin(1..$maxbin){
							#$hash{$md}[$group_number]{$bin}[0] = "NA";
							#$hash{$md}[$group_number]{$bin}[1] = "NA";
						#}
					#}
					#else {
					#	foreach my $bin(1..$maxbin){
					#		#print "in hash: ".$hash{$md}[$group_number]{$bin}[0]."\n";
					#		#print "norm ".$norms{$md}[$group_number]."\n";
					#		#print "sum ".$sums{$md}[$group_number]."\n";
					#		$hash{$md}[$group_number]{$bin}[0] = $hash{$md}[$group_number]{$bin}[0]*$norms{$md}[$group_number]/$sums{$md}[$group_number];
					#		$hash{$md}[$group_number]{$bin}[1] = $hash{$md}[$group_number]{$bin}[1]*$norms{$md}[$group_number]/$sums{$md}[$group_number];
					#	}
					#}
					
					#my $filehandle = $filehandles{$md}{$group_number};
					#print "going to print something\n";
					#foreach my $bin(1..$maxbin){
						#print $filehandle $hash{$md}[$group_number]{$bin}[0].",".$hash{$md}[$group_number]{$bin}[1].",";
					#}
					#print $filehandle "\n";
				
				
			#	}
			#}
			

			
			#$iteration_number++;
			#print $iteration_number."\n";
		}
		
		close GULP;
}

sub test_distr_to_stathist_norm {
	logic_global_median_statistics("h3", "nsyn", 0, "distrtest");
}

#test_obsv_manipulation();

sub test_obsv_manipulation {
	my $prot = "h1";
	my $restriction = 0;
	set_mutmap($prot, "nsyn");
	set_distance_matrix($prot);
	my %obs_vectors = get_observation_vectors();
	
	my %shuffled_obs_vectors = shuffle_observation_vectors(\%obs_vectors);
	my @mock_mutmaps = read_observation_vectors(\%shuffled_obs_vectors);
	%static_subs_on_node = %{$mock_mutmaps[0]};
	%static_nodes_with_sub = %{$mock_mutmaps[1]};

	my %hash = depth_groups_entrenchment_optimized(10,$restriction);
	
	
	open CSV, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_vector_simulation_test_".$restriction.".csv";
	foreach my $bin(1..32){
			print CSV $bin."_obs,".$bin."_exp,";
	}
	print CSV "\n";

	#foreach my $i(0..$iterations-1){
		foreach my $bin(1..32){
			print CSV $hash{$bin}[0].",".$hash{$bin}[1].",";
		}
		print CSV "\n";
	#}
	close CSV;
	
	
	#push @simulated_hists, \%hash;
	#%static_ring_hash = ();
	#%static_depth_hash = ();
	
	
}
sub test_remove_zero_branches {
	my $tree = parse_tree("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/h1.l.r.newick");
	my $multitree = remove_zero_branches($tree);
	open TREE, ">multitree.tre" or die "Cannot create file";
	my $treestr=tree2str($multitree);
	print TREE $treestr;
	close TREE;
}
