# syn research method; detecting ancestor syn mutations which can undergo several different syn muts
# and checking how many muts of each type really happened ({tv}GAG->GAT - 1, GAG->GAA = 0; {ts}GAG->GAA - 4)
sub reversals_list {
	my $self = shift;
	#my @group = @{$_[0]};
	my $root = $self->{static_tree}-> get_root;
	my @array;
	#unless (@group) {
	my @group = (1..$self->mylength());
	#}
		my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( 0.5, $root);
	$self->visitor_coat ($root, \@array,\&synresearch_visitor,\&no_check,\@args,0);
	my $filepath = File::Spec->catfile($self->{static_output_base}, $self->{static_protein}."_reversals_list");
	open my $file, ">$filepath" or die "Cannot open $filepath\n";
	print $file "ancestor,list,number";
	foreach my $ind (@group){
			foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node = ${$nod};
			if (!$node->is_terminal){

				if ($self->{static_subtree_info}{$node->get_name()}{$ind}{"list"}) {
						my $str = concat($ind, $node->get_name());
						my $counter;
						my $list;
						foreach my $nodename (@{$self->{static_subtree_info}{$node->get_name()}{$ind}{"list"}}){
							$counter++;
							$list = $list.";".$nodename;
						}
						$list = substr($list,1);
						$str = $str.",".$list.",".$counter."\n";
						print $file $str;
				}

			}
			
		}
	}
	$self->printFooter($file);
	close $file;
	
}







# syn research method; detecting ancestor syn mutations which can undergo several different syn muts
# and checking how many muts of each type really happened ({tv}GAG->GAT - 1, GAG->GAA = 0; {ts}GAG->GAA - 4)
sub synmut_types {
	my $self = shift;
	#my @group = @{$_[0]};
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my %results;
	#unless (@group) {
	my @group = (1..$self->mylength());
	#}
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( 0.5, $root);
	my $comparator = compare::new();
	$self->visitor_coat ($root, \@array,\&synresearch_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
			foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node = ${$nod};
			if (!$node->is_terminal){
				my $anc = $self->{static_subs_on_node}{$node->get_name()}{$ind}->{"Substitution::derived_allele"};
			#	print "anc $anc\n";
				my $synmuts = $comparator->get_synmuts($anc);
				my $numtypes = (scalar keys %{$synmuts->{"ts"}}) + (scalar keys %{$synmuts->{"tv"}});
				# todo create all possible syn types
				if ($numtypes > 1){
					if ($self->{static_subtree_info}{$node->get_name()}{$ind}{"synresearch"} && scalar @{$self->{static_subtree_info}{$node->get_name()}{$ind}{"synresearch"}} > 1) {
						my $str =  $ind."_".$node->get_name()." anc $anc\t";
						my $counter;
						foreach my $subst (@{$self->{static_subtree_info}{$node->get_name()}{$ind}{"synresearch"}}){
							if (defined $synmuts->{"ts"}{$subst->{"Substitution::derived_allele"}}){
								if ($synmuts->{"ts"}{$subst->{"Substitution::derived_allele"}} == 0 ) {$counter++;}
								$synmuts->{"ts"}{$subst->{"Substitution::derived_allele"}} +=1;
								$str = $str."ts sub ".$subst->{"Substitution::derived_allele"}."\t";
							}
							elsif (defined $synmuts->{"tv"}{$subst->{"Substitution::derived_allele"}}) {
								if ($synmuts->{"tv"}{$subst->{"Substitution::derived_allele"}} == 0 ) {$counter++;}
								$synmuts->{"tv"}{$subst->{"Substitution::derived_allele"}} +=1;
								$str = $str."tv sub ".$subst->{"Substitution::derived_allele"}."\t";
							}
						}
						$str = $counter."\t".$str."\n";
						print $str;
					}
				}
				$results{$ind}{$node->get_name()}{"anc"} = $anc;
				$results{$ind}{$node->get_name()}{"muts"} = $synmuts;
			}
			
		}
	}
	return %results;
	
}

  
  # entrenchment_visitorversion for for synmut_types: collect substitutions for each ancestor mutation (not just count them; do not track distance)
    	sub synresearch_visitor {
 		my $self = shift;
 		my $node = $_[0];
 		my $step = $_[1]->[0];
 		my $subtract_tallest = $self->{static_subtract_tallest};
 		my $no_neighbour_changing = $self->{static_no_neighbour_changing};
 		my $no_leaves = $self->{static_no_leaves};
 		my $comparator = $self->{static_comparator};
		if (!$node->is_root){
		my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")}; # closest_ancestors: ancestor mutation node for this node, key is a site number
		my $nname = $node->get_name();
		if (%closest_ancestors){
			foreach my $site_index(keys %closest_ancestors){ 
				my $anc_node = $closest_ancestors{$site_index};
				my $depth = get_sequential_distance($anc_node,$node);
				$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($depth,$step)}[1] += $node->get_branch_length;
			#	print "anc ".$anc_node->get_name()." site ".$site_index." node ".$node->get_name()." depth $depth bin ".bin($depth,$step)." branchlength ".$node->get_branch_length."\n";
				my $current_maxdepth = $self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth"};
				if ($current_maxdepth < $depth){
						$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth"} = $depth;
						if ($subtract_tallest){
							$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth_node"} = $nname;
						}
				}					
			}
			
			my @ancestors = keys %closest_ancestors;	
			foreach my $site_index(@ancestors){
				if (!($self->has_no_background_mutation($nname, $site_index))){
					delete $closest_ancestors{$site_index};
				}
			}	
		}
		
		foreach my $site_index(keys %{$self->{static_subs_on_node}{$nname}}){
			
			if ($closest_ancestors{$site_index}){
				my $anc_node = $closest_ancestors{$site_index};
				#my $halfdepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()} - ($node->get_branch_length)/2; #19.09.2016 
				#my $fulldepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()}; #19.09.2016 
			#	print " ancestor ".$anc_node->get_name(). " node ".$node->get_name()." depth $depth\n";
			#	push $static_subtree_info{$anc_node->get_name()}{$site_index}{"nodes"}, \$node;
				if (!$no_neighbour_changing || ($no_neighbour_changing && ! $comparator->is_neighbour_changing($self->{static_subs_on_node}{$nname}{$site_index}, 1))){
					if (!$no_leaves || ($no_leaves && !($node->is_terminal()))){
						push @{$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"synresearch"}}, $self->{static_subs_on_node}{$nname}{ $site_index}; 
						push @{$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"list"}}, $nname;
					}
				}
			#	$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($halfdepth,$step)}[1] += ($node->get_branch_length)/2; # #19.09.2016  15.09.2016 version: halves of branches with foreground mutations are trimmed (the only thing I changed here) 
			#	$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($fulldepth,$step)}[1] -= $node->get_branch_length; #19.09.2016 we added this length before, but shouldn't have done it
			#	print "mutation! anc ".$anc_node->get_name()." site ".$site_index." node ".$node->get_name()." depth $depth bin ".bin($depth,$step)." branchlength ".$node->get_branch_length."\n";
				
			}
			$closest_ancestors{$site_index} = $node;
			
		}
		
		$node->set_generic("-closest_ancestors" => \%closest_ancestors);
		}
#!
 	} 