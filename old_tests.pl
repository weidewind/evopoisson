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