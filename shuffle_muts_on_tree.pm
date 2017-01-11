#This module provades methods for shuffling mutations on the tree according to Poisson like process
use strict;
use Bio::Phylo::IO;
use List::BinarySearch qw(binsearch_pos);
use Class::Struct;
use Math::Random::ISAAC;

my $rng = Math::Random::ISAAC->new(localtime());
my $MaxTries=1;

struct Strip =>{
	nodes => '@',
	square => '$'
}

struct StripConstrains =>{
	number_of_mutations => '$',
	lifetime => '$',
	stoppers => '@'
}

#This function partitions tree into nonoverlapping strips
sub strip_tree{
	my ($rnode,$ra_time_stamps,$ra_out_strips)=@_;
	my $rtime=$rnode->get_generic('time');
	die "\nError strip_tree(): time stamps on nodes are expected!" unless defined $rtime;
	@{$ra_out_strips}=();
	$rnode->visit_breadth_first(
		-in   => sub{
			my $node=shift;
			my $time=$node->get_generic('time');
			die "\nError strip_tree(): time stamps on nodes are expected!" unless defined $time;
			$time-=$rtime;
			my $idx=binsearch_pos {$a<=>$b} $time,@{$ra_time_stamps};
			if($idx<@{$ra_time_stamps}){
				$ra_out_strips->[$idx]=Strip->new() unless defined $ra_out_strips->[$idx];
				push @{$ra_out_strips->[$idx]->nodes()},$node;
				my $sqr=$ra_out_strips->[$idx]->square();
				$ra_out_strips->[$idx]->square($sqr+$node->get_branch_length());
			}
		}
	);	
}

#The function expects a lifetime constrain in terms of a number of strips
#Number of strips is required!
sub shuffle_mutations{
	my ($rnode,$ra_strips,$ra_strip_constr,$ra_out_event_nodes)=@_;
	die "\nError shuffle_mutations(): all arrays with constrains are required to be of equal size!" unless scalar(@{$ra_mut_number})==scalar(@{$ra_strip_number});
	#stores detailed information about branch to strip assignment
	my %node2strip;
	#make cdf for distributing of mutations across strips
	my @CDF;
	$CDF[0]=0;
	for(my $i=0;$i<@{$ra_strips};$i++){
		$CDF[$i]=$ra_strips->[$i]->square();
		$CDF[$i]+=$CDF[$i-1] if $i;
		for(my $j=0;$j<@{$ra_strips->[$i]->nodes()};$j++){
			my $node=$ra_strips->[$i]->nodes($j);
			$node2strip{$node->get_name()}=[($i,$j)];
		}
	}
	@{$ra_out_event_nodes}=();
	push @{$ra_out_event_nodes},[];
	my $ii=$MaxTries;
	my $I=0;
	while($I<@{$ra_strip_constr}){
		my $n=$ra_strip_constr->[$I]->number_of_mutations();
		my $strip_number=$ra_strip_constr->[$I]->lifetime(); #number of strips
		die "\nError shuffle_mutations(): wrong number of strips!" unless $strip_number>0&&$strip_number<=@{$ra_strips};
		my $ra_events=$ra_out_event_nodes->[-1];
		my @nsamples;
		my %blocked;
		if(defined $ra_strip_constr->[$I]->stoppers()){
			foreach my $node(@{$ra_strip_constr->[$I]->stoppers}){
				$node->visit_breadth_first(
					-in   => sub{
						my $nd=shift;
						$blocked{$nd->get_name}=1;
					}
				);
			}
		}
		my %mutations;
		#seed mutations across strips
		while($n--){
			$_=$rng->rand()*$CDF[$strip_number-1];
			my $idx=binsearch_pos {$a<=>$b} $_,@CDF;
			$nsamples[$idx]++;
		}
		#look through strips starting from the furthest
		my $i=$#nsamples;
		for(;$i>-1;$i--){
			#place mutations on tree branches within a strip
			next unless defined $nsamples[$i];
			my $n=$nsamples[$i];
			#make cdf for the current strip
			my @cdf;
			$cdf[0]=0;
			for(my $j=0;$j<@{$ra_strips->[$i]->nodes()};$j++){
				my $node=$ra_strips->[$i]->nodes($j);
				if(!exists $blocked{$node->get_name()}){
					$cdf[$j]=$node->get_branch_length();
				}
				$cdf[$j]=$cdf[$j-1] if $j;
			}
			while($n){
				last unless $cdf[0]<$cdf[-1];
				#sample branch in the current strip for the next mutation
				$_=$rng->rand()*$cdf[-1];
				my $idx=binsearch_pos {$a<=>$b} $_,@cdf;
				my $node=$ra_strips->[$i]->nodes($idx);
				my $pnode=$node;
				#check if the sampled brannch is plausible
				while($pnode!=$rnode){
					my $pname=$pnode->get_name;
					if(exists $blocked{$pname}){
						$pnode=$rnode unless exists $mutations{$pname};
						last;
					}else{
						$blocked{$pname}=1;
						#recalculate @cdf
						my ($pstrip_idx,$pidx)=($node2strip{$pname}->[0],$node2strip{$pname}->[1]);
						if($pstrip_idx==$i){ #$pnode is in the current strip 
							my $l=$ra_strips->[$i]->nodes($pidx)->get_branch_length(); 
							for(my $j=$pidx;$j<@{$ra_strips->[$i]->nodes()};$j++){
								$cdf[$j]-=$l;
							}
						}
					}
					$pnode=$node->get_parent;
				}
				if($pnode==$rnode){
					#mutation could be placed on the sampled branch
					$mutations{$node->get_name}=1;
					push @{$ra_events},$node;
					$n--;
				}#else try again
			}
			if($n){
				#failed to place mutations in the current strip
				last;
			}
		}
		if($i>-1){
			#failed to place mutatons
			@{$ra_events}=();
			if($ii--==0){ 
				#no more attempts allowed
				$ra_out_event_nodes->[-1]=undef;
				$i=-1;
			}
		}
		if($i==-1){
			#go next
			push @{$ra_out_event_nodes},[];
			$ii=$MaxTries;
			$I++;
		}
	}
}

#input args: 
#0) Bio::Phylo::IO tree;
#1) hash{node_name}{site_index} = StripConstrains object
#
#output args:
#hash{node_name}{site_index} = [], ref on array with nodes carrying mutations 
sub shuffle_mutations_on_tree{
	my ($tree_,$rh_constrains,$rh_out_subtree)=@_;
	%{$rh_out_subtree}=();
	#if a setting of attributes on tree nodes is not a problem than cloning may be omitted
	my $tree=$tree_->clone();
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
	#making time stamps to partition branches of corresponding subtrees of nodes into strips
	#width of strips are defined by the dencity of sampling of strains
	my %timestamps;
	foreach my $node(@{$tree->get_internals}){
		my $name=$node->get_name;
		$timestamps{$name}=[];
	}
	$tree->visit_depth_first(
		-in   => sub{
			my $node=shift;
			my $name=$node->get_name;
			my $time=$node->get_generics('time');
			if(!$node->is_root){
				my $pnode=$node->get_parent;
				my $pname=$pnode->get_name
				if($node->is_terminal){
					push @{$timestamps{$pname}},$time;
				}else{
					push @{$timestamps{$pname}},@{$timestamps{$name}};
				}
			}
			my @srt=sort {$a<=>$b} @{$timestamps{$name}};
			@{$timestamps{$name}}=();
			push @{$timestamps{$name}},$srt[0]-$time;
			for(my $i=1;$i<@srt;$i++){
				push @{$timestamps{$name}},$srt[$i]-$time if $srt[$i]>$str[$i-1];
			}	
		}
	);
	foreach my $node(@{$tree->get_internals}){
		my $name=$node->get_name;
		if(defined $rh_constrains->{$name}){
			$rh_out_subtree->{$name}={};
			my $max_life_time=0;
			my @ts=@{$timestamps{$name}};
			#trunkate time stamp array by the maximal life time over all sites on the branch
			my $n=keys %{$rh_constrains->{$name}}; #number of sites on a branch
			foreach my $site(keys %{$rh_constrains->{$name}}){
				my $strc=$rh_constrains->{$name}->{$site};
				if(defined $strc->lifetime()){
					my $l=$strc->lifetime();
					$max_life_time=$l if $max_life_time<$l;
					$n--;
				}
			}
			$max_life_time=0 if $n; 
			if($max_life_time){
				my $idx=binsearch_pos {$a<=>$b} $max_life_time,@ts;
				$#ts=$idx unless $idx>$#ts;
			}
			my @strips;
			strip_tree($node,\@ts,\@strips);
			my @strip_constrs;
			my @sites;
			foreach my $site(keys %{$rh_constrains->{$name}}){
				push @sites,$site;
				my $strc=StripConstrains->new();
				$strc->number_of_mutations($rh_constrains->{$name}->{$site}->number_of_mutations);
				$strc->stoppers($rh_constrains->{$name}->{$site}->stoppers);
				#express a lifetime as a number of strips
				$max_life_time=0;
				if(defined $rh_constrains->{$name}->{$site}->lifetime()){
					$max_life_time=$rh_constrains->{$name}->{$site}->lifetime();
				}
				my $n=@ts; #maximal number of strips
				if($max_life_time){
					$n=binsearch_pos {$a<=>$b} $max_life_time,@ts;
				}
				$strc->lifetime($n);
				push @strip_constrs,$strc;
			};
			my @tmp; 
			shuffle_mutations($node,\@strips,\@strip_constrs,\@tmp);
			for(my $i=0;$i<@sites;$i++){
				my $site=$sites[$i];
				if(defined $tmp[$i]){
					$rh_out_subtree->{$name}->{$site}=[] unless defined $rh_out_subtree->{$name}->{$site};
					die "\nNot all mutations were positioned on the tree!" unless $mnumber[$i]==scalar(@{$tmp[$i]});
					push @{$rh_out_subtree->{$site}},@{$tmp[$i]};
				}else{
					$rh_out_subtree->{$name}->{$site}=undef;
				}
			}
		}		
	}
}

