#This module provades methods for shuffling mutations on the tree according to Poisson like process
package shuffle_muts_on_tree_exp;

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use Math::Random::ISAAC;

my $rng = Math::Random::ISAAC->new(localtime());
my $MaxTries=100;

struct Constrains =>{
	number_of_mutations => '$',
	hazard => '$',
	stoppers => '@'
};


#The function expects a lifetime constrain in terms of a number of strips
#Number of strips is required!
sub shuffle_mutations {
	my ($rnode,$ra_strip_constr,$ra_out_event_nodes, $poisson)=@_;
	@{$ra_out_event_nodes}=();
	push @{$ra_out_event_nodes},[];
	my $ii=$MaxTries;
	my $I=0;
	while($I<@{$ra_strip_constr}){
		my $n=$ra_strip_constr->[$I]->number_of_mutations();
		my $lambda=$ra_strip_constr->[$I]->hazard(); #number of strips
		print "-----------------------ancestor is ".$rnode->get_name."site $I, want $n mutations here \n";
		print " lambda for ".$rnode->get_name." anc site $I is $lambda\n";
		my $ra_events=$ra_out_event_nodes->[-1];
		my @nsamples;
		my %blocked;
		#if(defined $ra_strip_constr->[$I]->stoppers()){
		#	foreach my $node(@{$ra_strip_constr->[$I]->stoppers}){
		#		print "Stopper here ".$node->get_name." \n";
		#		$blocked{$node->get_name}=1;
		#	}
		#}
		if(defined $ra_strip_constr->[$I]->stoppers()){
			NODE: foreach my $node(@{$ra_strip_constr->[$I]->stoppers}){ 
				next if ($node->get_name() eq $rnode->get_name()); # ignore stopper if it is on the same branch (may happen in fake samples)
				$blocked{$node->get_name()}=1;
				#print " already blocked for ".$rnode->get_name()." : ".$node->get_name()."\n";
				if (!($node ->is_terminal)) { #added
					my $child = $node->get_child(0); # added: visit_ also visits starting node's sister
					$child->visit_breadth_first(
						-in   => sub{
							my $nd=shift;
							next NODE if ($nd->get_name() eq $rnode->get_name()); # ignore stopper if it is higher than the ancestor branch
							$blocked{$nd->get_name()}=1;
							#print " already blocked for ".$rnode->get_name()." : ".$nd->get_name()."\n";
						}
					);
				}
			}
		}

		#seed mutations following exponential distribution
		my $head=0;
		my @stack;
		push @stack,$rnode unless exists $blocked{$rnode->get_name};
		my $t0=$rnode->get_generic('time');
		#while($head<@stack){
		while(@stack){
			#my $node=$stack[$head++];
			my $node= pop @stack;
			if($node!=$rnode && !$blocked{$node->get_name}){
				my $pnode=$node->get_parent;
				#my $t=$pnode->get_generic('time')-$t0;
				#my $p=exp(-$lambda*$t);
				#$t=$node->get_generic('time')-$t0;
				#$p-=exp(-$lambda*$t);
				my $t=$node->get_generic('time')-$pnode->get_generic('time');
				my $p;
				if($poisson){
					$p=$lambda*$t;
				}
				else {
					$p=1-exp(-$lambda*$t);
					print "for ".$node->get_name." time is ".$t." and p is $p \n";
				}
				$_=$rng->rand();
				if($_<=$p){ 
					push @{$ra_events},$node->get_name;
					#print "chose node ".$node->get_name." for anc ".$rnode->get_name."\n";
					#print "Now get ";
					#foreach my $nname(@{$ra_events}){
					#	print $nname."\t";
					#}
					#print "\n";
					$blocked{$node->get_name}=1;
					$n--;
				}
			}
			if(!$blocked{$node->get_name}){ 
				#print " going to push ".$node->get_name." children to stack\n";
				foreach my $chnode(@{$node->get_children}){ 
					push @stack,$chnode unless exists $blocked{$chnode->get_name};
					#print "pushed ".$chnode->get_name." to stack \n" unless exists $blocked{$chnode->get_name};
				}
			}
		}
		if($n){
			#print "n is $n \n";
			#failed to place all mutations
			@{$ra_events}=();
			if($ii--==0){ 
				#no more attempts allowed
				print " failed to place ".$ra_strip_constr->[$I]->number_of_mutations()." mutations at ".$rnode->get_parent->get_name."\n";
				$ra_out_event_nodes->[-1]=undef;
			}
		}
		if($n==0 ||$ii==-1){ # changed: || $ii==0
			#go next
		#	print "lucky me! \n";
			push @{$ra_out_event_nodes},[];
			$ii=$MaxTries;
			$I++;
		}
	}
}

#input args: 
#0) Bio::Phylo::IO tree;
#1) hash{node_name}{site_index} = Constrains object
#
#output args:
#hash{node_name}{site_index} = [], ref on array with nodes carrying mutations 
sub shuffle_mutations_on_tree{
	my ($tree,$rh_constrains, $poisson)=@_;
	my $rh_out_subtree;
#	print "Went inside shuffler\n";
	#if a setting of attributes on tree nodes is not a problem than cloning may be omitted
	#my $tree=$tree_->clone();
	#print "Tree cloned\n";
	#set time on the tree nodes
	$tree->visit_breadth_first(
		-in   => sub{
			my $node=shift;
			if($node->is_root){
				$node->set_generic('time' => 0);
			}else{
				my $pnode=$node->get_parent;
				my $time=$pnode->get_generic('time');
				$time+=$node->get_branch_length;
				$node->set_generic('time' => $time);
			}
		}
	);	
	#print " Tree visited\n";
	#my @sites;
	#my @strip_constrs;
	foreach my $node(@{$tree->get_internals}){
		my @sites;
		my @strip_constrs;
		my $name=$node->get_name;
		if(defined $rh_constrains->{$name}){
			$rh_out_subtree->{$name}={};
			foreach my $site(keys %{$rh_constrains->{$name}}){
				push @sites,$site;
				print "pushed $site into site array for $name\n";
				push @strip_constrs,$rh_constrains->{$name}->{$site};
			};
			my @tmp; 
			print " just about to shuffle muts on tree for node $name\n";
			shuffle_mutations($node,\@strip_constrs,\@tmp, $poisson);
			for(my $i=0;$i<@sites;$i++){
				my $site=$sites[$i];
				if(defined $tmp[$i]){
					$rh_out_subtree->{$name}->{$site}=[] unless defined $rh_out_subtree->{$name}->{$site};
					print "for $site $name need ".$strip_constrs[$i]->number_of_mutations." mutations, got ".scalar(@{$tmp[$i]})."\n";
					die "\nNot all mutations were positioned on the tree!" unless $strip_constrs[$i]->number_of_mutations==scalar(@{$tmp[$i]});
					push @{$rh_out_subtree->{$name}->{$site}},@{$tmp[$i]};
				}else{
					$rh_out_subtree->{$name}->{$site}=undef;
				}
			}
		}		
	}
	return $rh_out_subtree;
}

1;