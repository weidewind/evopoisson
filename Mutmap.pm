#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

package MutMap;

use strict;
use Bio::Phylo::IO;
use DnaUtilities::compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing);
use TreeUtils::Phylo::PhyloUtils qw(remove_zero_branches);

use Bio::Tools::CodonTable;
use TreeUtils::Phylo::FigTree;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Try::Tiny;
use List::Util qw(sum);
use Const::Fast;
	use Switch;

use List::Util qw/shuffle/; 
use Statistics::Basic qw(:all);
use Statistics::TTest;
use Statistics::Descriptive;
use Storable qw(store retrieve lock_retrieve);

use Class::Struct;
use DnaUtilities::observation_vector qw(make_observation_vector shuffle_obsv);
use IO::Handle;
use Cwd qw(abs_path cwd getcwd);
use File::Spec;
use Clone 'clone';
use Sub::Identify ':all';
use autodie;

$| = 1;

	
	
#rewrite
	sub set_node_distance {
		$distance_hash{$_[0]}->{$_[1]} = $_[2];	
	}
#rewrite	
	sub set_alignment_length {
		$static_alignment_length = $_[0]; 
	}
#rewrite	
	sub has_node_distance {
		if (!defined $distance_hash{$_[0]}->{$_[1]}){
			return 0;
		}
		else {
			return 1;
		}
	}
#rewrite	
	sub get_node_distance {
		return $distance_hash{$_[0]}->{$_[1]};
	}
	
	
	struct Mutation => {
		site_index => '$',
		node => '$',
	};
	

	
	sub pathFinder {
		my $mutmap = shift;
		my $subroutine = shift;
		my $bigtag = shift;
		my $bigdatatag = shift;
		
		#my $sub_name = (caller(0))[3];
		my $subroutine_name = split(\::\, $subroutine)[-1];
		
		my $output_base = File::Spec->catdir(getcwd(), "output", $bigdatatag, $bigtag, state_tag($mutmap -> {state}), maxpath_tag($mutmap -> {subtract_tallest}), $mutmap -> {protein}); 
		my $input_base = File::Spec->catdir(getcwd(), "data", $bigdatatag,);
	}
	
	sub maxpath_tag{
		my $subtract_maxpath = $_[0];
		my $tag;
		if ($subtract_maxpath){
			$tag = "maxpath_subtracted";
		}
		else {
			$tag = "maxpath_not_subtracted"
		};
		return $tag;
	}


	sub state_tag {
		my $state = $_[0];
		my $tag;
		if ($state eq "s" || $state eq "syn") { $tag = "syn";}
		elsif ($state eq "n" || $state eq "nsyn") { $tag = "nsyn";}
		else {die "Unknown state $state; expected syn or nsyn";}
		return $tag;
	}
	
	sub syn_tag {
		my $syn = $_[0];
		my $tag;
		if ($syn == 1){
			$tag = "syn";
		}
		else {
			$tag = "nsyn";
		}
		return $tag;
	}

	
	
	sub new {
		my ($class, $args) = @_;	
		my $output_base = File::Spec->catdir(getcwd(), "output", $args->{bigdatatag}, $args->{bigtag}, state_tag($args->{state}), maxpath_tag($args->{subtract_tallest}), $args->{protein}); 
		my $input_base = File::Spec->catdir(getcwd(), "data", $args->{bigdatatag},);
		my $treefile = File::Spec->catfile($input_base, $args->{protein}.".l.r.newick");
		my $static_tree = parse_tree($treefile)  or die "No tree at $treefile";
		
		if ($args->{fromfile}){
			my $realdatapath = File::Spec->catfile($output_base, $args->{protein}."_realdata");
			my $realdata = retrieve ($realdatapath) or die "Cannot retrieve ".$realdatapath;
			
			my $self = { 
				static_output_base => $output_base,
				static_input_base => $input_base,
				static_protein => $args->{protein},
				static_subtract_tallest => $args->{subtract_tallest},
				static_tree => $static_tree,
				static_state => $args->{state},
				static_hash_of_nodes => $realdata->{"hash_of_nodes"}, 
				static_distance_hash => $realdata->{"distance_hash"},
				static_subs_on_node => $realdata->{"subs_on_node"}, # we never use these two when we produce new mutmappers from file (they are taken from observaton_vectors)
				static_nodes_with_sub => $realdata->{"static_nodes_with_sub"}, #
				static_background_subs_on_node => $realdata->{"bkg_subs_on_node"},
				static_background_nodes_with_sub => $realdata->{"bkg_nodes_with_sub"},
				realdata => $realdata,
			};
		}
		else {
			my @arr = parse_fasta(File::Spec->catfile($input_base, $args->{protein}.".all.fa"));
			my %fasta = %{$arr[0]};
			my $alignment_length = $arr[1];
			my $static_protein  = $args->{protein};
			my %static_fasta = %fasta;
			my %static_hash_of_nodes;	
			my @nodes = $static_tree -> get_nodes;
			foreach my $node(@nodes){
				#if ($node->is_root()) {next;}
				my $name = $node ->get_name();
				$self ->{static_hash_of_nodes}{$name} = \$node;
			}
			
			my @mutmaps;
			my @bkg_mutmaps;
			if($args->{state} eq "syn"){
				@mutmaps = synmutmap($static_tree, \%fasta);
				@bkg_mutmaps = codonmutmap($static_tree, \%fasta);
			} 
			elsif($args->{state} eq "nsyn"){
				@mutmaps = codonmutmap($static_tree, \%fasta);
				@bkg_mutmaps = synmutmap($static_tree, \%fasta);
			} 
			else {
				die "only syn or nsyn can be used as the second argument; unknown ."$args->{state}." was used instead";
			}

			my $self = {
				static_output_base => $output_base,
				static_input_base => $input_base,
				static_protein => $args->{protein},
				static_alignment_length => $alignment_length, 
				static_subtract_tallest => $args->{subtract_tallest},
				static_tree => $static_tree,
				static_fasta => { %static_fasta },
				static_state  => $args->{state},
				static_hash_of_nodes => { %static_hash_of_nodes },
				static_subs_on_node => $mutmaps[0],
				static_nodes_with_sub => $mutmaps[1],
				static_background_subs_on_node => $bkg_mutmaps[0],
				static_background_nodes_with_sub => $bkg_mutmaps[1],
			};
		}	

		bless $self, $class;
		return $self;
	}
	
	





## according to http://www.biomedcentral.com/1471-2148/10/253
		const my @n1_decreasing => ("AGG", "TCG", "GAT", "CGT", "ACC", "GCC", "CAG", "GGG", "GGC");

		const my @h1_decreasing => ("ACG", "TCA", "CTC", "GCG", "GCA", "CCG", "TGC", "GTG");

		const my @n2_decreasing => ("AAT", "CTC", "GAG", "TCT", "ACT", "TGT", "CCG", "GGC", "GAC", "AAA", "TCA");

		const my @h3_decreasing => ("CTG", "CGC", "CCT", "TGC",  "GAC", "AGG", "TAT", "AAG", "GGG", "CGG");

		const my @n1_increasing => ("AGA", "ACA", "GGA", "CAC", "TCA", "CTT", "CAA", "AGT");

		const my @h1_increasing  => ("ACA", "GCC", "CCT", "TGT", "TCC", "AGC");

		const my @n2_increasing => ("AAC", "TCC", "GAA", "GTT", "TGC", "GAT", "AAG", "GCC", "ACA");
		
		const my @h3_increasing => ("TTG", "AGA", "TGT", "GCC",  "CTA", "GAT", "TAC", "CCG", "GGA", "AAA", "CCC");

		const my @all_codons => ("TCA", "TCC", "TCG", "TCT", "TTC", "TTT", "TTA", "TTG", "TAC", "TAT", "TAA", "TAG", "TGC", "TGT", "TGA", "TGG", "CTA", "CTC", "CTG", "CTT", "CCA", "CAT", "CAA", "CAG", "CGA", "CGC", "CGG", "CGT", "ATA", "ATC", "ATT", "ATG", "ACA", "ACC", "ACG", "ACT", "AAC", "AAT", "AAA", "AAG", "AGC", "AGT", "AGA", "AGG", "CCC", "CCG", "CCT", "CAC", "GTA", "GTC", "GTG", "GTT", "GCA", "GCC", "GCG", "GCT", "GAC", "GAT", "GAA", "GAG", "GGA", "GGC", "GGG", "GGT");

		

## returns a hash: key - codon, value - -1, if it is decreasing over time,  
##                                       1, if it is increasing, 
##										 0 otherwise.
sub codon_evolution{
	my $protein = $_[0];
	my %hash;
	@hash{@all_codons} = 0;
	switch($protein){
		case "n1" {	@hash{@n1_decreasing} = -1; 
			        @hash{@n1_increasing} = 1; }
		case "n2" {	@hash{@n2_decreasing} = -1; 
					@hash{@n2_increasing} = 1; }
		case "h1" {	@hash{@h1_decreasing} = -1; 
					@hash{@h1_increasing} = 1; }
		case "h3" {	@hash{@h3_decreasing} = -1; 
					@hash{@h3_increasing} = 1; }
		}
	
	return %hash;
}



# read .fasta into a hash
sub parse_fasta {
	my $nodeseqs_file = shift;
	my %nodeseqs;
	my $seqio = Bio::SeqIO->new(-file => $nodeseqs_file, -format => "fasta");
	my $length;
	while ( my $seqobj = $seqio->next_seq ) {
		my $trimmed_id = (split(/\//, $seqobj->display_id))[0];
    	$nodeseqs{ $trimmed_id } = $seqobj->seq;
    	if (!$length){
    		$length = $seqobj->length();
    	}
	}
	return (\%nodeseqs, $length);
}

sub parse_tree {
					my $tree_file = $_[0];
					open TREE, "< $tree_file" or die "Cannot open file ".$tree_file."\n";
					# get a newick string from some source
 					my $tree_string;
 					my $t;
 					while($t = <TREE>){
 						$t =~ s/\n$//;
 						$tree_string = $tree_string.$t;
 					}
 					 close TREE;

 					# Call class method parse from Bio::Phylo::IO
 					# note: newick parser returns 'Bio::Phylo::Forest'
                    # Call ->first to retrieve the first tree of the forest.
 					my $tree = Bio::Phylo::IO->parse(
  					  -string => $tree_string,
   					  -format => 'newick'
 					)->first;

 					return $tree;
	
} 


sub codonmutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		my %nsyn = nsyn_substitutions_codons($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
};


sub synmutmap {

	my $tree = $_[0];
	my %nodeseqs = %{$_[1]};
	my %subs_on_node;
	my %nodes_with_sub;
	
	my @nodes = $tree -> get_nodes;
	foreach my $node(@nodes){
		if ($node->is_root()) {next;}
		my $name = $node ->get_name();
		my %nsyn = syn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
									  $nodeseqs{$name});					  
		$subs_on_node{$name}=\%nsyn;
		for	my $site_index(keys %nsyn){
			if (! exists $nodes_with_sub{$site_index}){
				$nodes_with_sub{$site_index} = ();
			}
			push (@{$nodes_with_sub{$site_index}}, \$node);
		}	
	}

	return (\%subs_on_node, \%nodes_with_sub);
}


# sets static_sorted_nodnames and static_sorted_sites, retruns incidence_hash
sub incidence_matrix {
	my $self = shift;	
	my %matrix;
	my $length;
	if ($self->{static_state} eq "nsyn"){
		$length = $self->{static_alignment_length}/3;
	}
	elsif ($self->{static_state} eq "syn") {
		$length = $self->{static_alignment_length};
	}
	
	my @sorted_sites;
	my @sorted_nodnames;
	
	# select nodes with at least one mutation of the corresponding type (syn or nsyn, depending on the mutmap state)
	my @nodes = $self -> {static_tree} -> get_nodes;
	foreach my $node(@nodes){
		my $name = $node ->get_name();
		if (scalar (keys %{$self->{static_subs_on_node}{$name}}) > 0){
			push @sorted_nodnames, $name;
		}
	}
	
	my %empty_nodes_hash = map { $_ => 0 } @sorted_nodnames;
	my %incidence_hash;
	
	# select sites with at least 3 mutations of the corresponding type
	# upd - not sure if such sites should be excluded, changed minimum to 1
	foreach my $ind(1..$length){
		if($self->{nodes_with_sub}{$ind} && scalar @{$self->{nodes_with_sub}{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
			my %site_incidence = %empty_nodes_hash;
			push @sorted_sites, $ind;
			#print " added $ind to sorted sites\n";
			foreach my $node(@{$self->{nodes_with_sub}{$ind}}){
				$site_incidence{$$node->get_name()} = 1;
			}
			$incidence_hash{$ind} = \%site_incidence;
		}
	}
	
	$self->{static_sorted_nodnames} = \@sorted_nodnames;
	$self->{static_sorted_sites} = \@sorted_sites;
	
	return %incidence_hash;
};



sub print_incidence_matrix {
	my $self = $_[0];
	my %incidence_hash = %{$_[1]};
	#my $path = $_[2];
	my $statetag = state_tag($self->{static_state});
	my $matrix_file = File::Spec -> catfile($self->{static_output_base}, $self->{static_protein}."_".$statetag."_incidence_matrix");
	my $sorted_sites_file = File::Spec -> catfile($self->{static_output_base}, $self->{static_protein}."_".$statetag."_sorted_sites");
	my $sorted_nodnames_file = File::Spec -> catfile($self->{static_output_base}, $self->{static_protein}."_".$statetag."_sorted_nodnames");
	
	# todo check if static_sorted_nodnames exists, if not - throw error and die (where did you take incidence_hash from?)
	unless ($self->{static_sorted_nodnames} && $self->{static_sorted_sites}) die "There is no static_sorted_nodnames (or static_sorted_sites) in this mutmapper. Where did you take incidence_hash from?\n";
	
	open MATRIX, ">$matrix_file" or die "Cannot open file ".$matrix_file."\n";
	foreach my $nodname (@{$self->{static_sorted_nodnames}}){
		foreach my $ind (@{$self->{static_sorted_sites}}){
			print MATRIX $incidence_hash{$ind}->{$nodname};
		}
		print MATRIX "\n";
	}
	close MATRIX;
	
	open SSITES, ">$sorted_sites_file" or die "Cannot open file ".$sorted_sites_file."\n";
	foreach my $ind(@{$self->{static_sorted_sites}}){
		print SSITES $ind."\n";
	}
	close SSITES;
	
	open SNODES, ">$sorted_nodnames_file" or die "Cannot open file ".$sorted_nodnames_file."\n";
	foreach my $name(@{$self->{static_sorted_nodnames}}){
		print SNODES $name."\n";
	}
	close SNODES;
	
}


sub read_incidence_matrix {
	my $self = $_[0];
	my $matrix_file = $_[1];
	open MATRIX, "<$matrix_file" or die "Cannot open file ".$matrix_file."\n";
	my %subs_on_node;
	my %nodes_with_sub;
	my $line_index = 0;
	while(<MATRIX>){
		if (/^$/) {last;}
			my $nodname = $self ->{static_sorted_nodnames}[$line_index];
			my @sites = split(',');
			my %substs;
			foreach my $s(@sites){
				my $ind = $self ->{static_sorted_sites}[$s-1];
#				print " $s is $ind\n";
				my $p=Substitution->new();
				$p->position($ind);
				$p->ancestral_allele("ATG");
				$p->derived_allele("ATG");
				$substs{$ind} = $p;
				if (! exists $nodes_with_sub{$ind}){
					$nodes_with_sub{$ind} = ();
				}
	#			print "\n nodndame ".$_."\n";
	#			print "\nREF 1 ".ref($static_hash_of_nodes{$nodname})."\n";
	#			print "\nREF 2 ".ref(${$static_hash_of_nodes{$nodname}})."\n";
				push (@{$nodes_with_sub{$ind}}, \${$self ->{static_hash_of_nodes}{$nodname}}); #вытащить из дерева по имени
			}
			$subs_on_node{$nodname} = \%substs;
			$line_index++;

	}
	close MATRIX;
	return (\%subs_on_node, \%nodes_with_sub);
			
}


## unlike original bio::phylo get_mrca, returns the node n1 closest to the root, if you give it two sequential nodes n1 and n2
## (in that case get_mrca returns the youngest ancestor of n1)

sub get_mrcn {
        my ( $node, $other_node ) = @_;
        if ( $node->get_id == $other_node->get_id ) {
            return $node;
        }
        my $self_anc  = $node->get_ancestors;
		unshift @{$self_anc}, $node;
        my $other_anc = $other_node->get_ancestors;
		unshift @{$other_anc}, $other_node;
        for my $i ( 0 .. $#{$self_anc} ) {
            my $self_anc_id = $self_anc->[$i]->get_id;
            for my $j ( 0 .. $#{$other_anc} ) {
                if ( $self_anc_id == $other_anc->[$j]->get_id ) {
                    return $self_anc->[$i];
                }
            }
        }

        return $self_anc->[-1];
    }

## the only difference from calc_patristic_distance is that it uses get_mrcn instead of get_mrca
## If you give it two sequential nodes, it returns the distance between them 

 sub calc_true_patristic_distance {
        my ( $node, $other_node ) = @_;
        my $patristic_distance = 0;
        my $mrca    = get_mrcn($node, $other_node);
        my $mrca_id = $mrca->get_id;
        while ( $node->get_id != $mrca_id ) {
            my $branch_length = $node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $node = $node->get_parent;
        }
        while ( $other_node and $other_node->get_id != $mrca_id ) {
            my $branch_length = $other_node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $other_node = $other_node->get_parent;
        }
        return $patristic_distance;
    }
    
    
## like calc_patristic_distance, but returns 0 for two sequential nodes    
    sub calc_my_distance {
        my ( $node, $other_node ) = @_;
        my $patristic_distance = 0;
        my $mrca    = get_mrcn($node, $other_node);
        my $mrca_id = $mrca->get_id;
        if ( $node->get_id == $mrca_id || $other_node->get_id == $mrca_id){
        	return 0;
        }
        while ( $node->get_id != $mrca_id ) {
            my $branch_length = $node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $node = $node->get_parent;
        }
        while ( $other_node and $other_node->get_id != $mrca_id ) {
            my $branch_length = $other_node->get_branch_length;
            if ( defined $branch_length ) {
                $patristic_distance += $branch_length;
            }
            $other_node = $other_node->get_parent;
        }
        return $patristic_distance;
    }
    
    sub node_distance {
    	my ( $node, $other_node ) = @_;
    	if  (has_node_distance($node, $other_node)){
    		return get_node_distance($node, $other_node);
    	}
    	else {
    		## calc_true instead of calc_my since 02 06 2015
    		my $dist = calc_true_patristic_distance($node, $other_node);
    		set_node_distance($node, $other_node, $dist);
    		return $dist;
    	}
    	
    }







# prints nexus tree, on wich all mutations in the specified site are shown 

sub print_static_tree_with_mutations{
	my $self = shift;
	my $site = shift;
	my $myCodonTable   = Bio::Tools::CodonTable->new();

	my %sites;
	my %color;
	foreach my $n(@{$self -> {static_nodes_with_sub}{$site}}){

		my $sub = ${$self -> {static_subs_on_node}{$$n->get_name()}}{$site};
		$sites{$$n->get_name()} = $$n->get_name()."_".$site."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
		$color{$$n->get_name()} = "-16776961";
	}

	open TREE, ">".$self -> {static_protein}."_sites_".$site.".tre";
	print TREE "#NEXUS\n\nbegin trees;\n";
	print TREE "\ttree $site = [&R] ";
	my $tree_name=tree2str($self -> {static_tree},sites => \%sites, color=>\%color);
	print TREE $tree_name;
	print TREE "\nend;";
	close TREE;

	foreach my $trn(@{$self -> {static_nodes_with_sub}{$site}}){
		print ($$trn->get_name()."\t");
		print (${$self -> {static_subs_on_node}{$$trn->get_name()}}{$site}->{"Substitution::derived_allele"});
		print "\n";
		foreach my $trr(@{$self -> {static_nodes_with_sub}{$site}}){
		
			print "\t".calc_true_patristic_distance($$trr, $$trn)."_";
			print (${$self -> {static_subs_on_node}{$$trr->get_name()}}{$site}->{"Substitution::derived_allele"}."_");
			print $$trr->get_name();
			print "\n";
		}
		print "\n";
	}
}




sub mean_ignore_nulls{
	if (!defined $_[0]){
		return 0;
	}

	my @arr = @{$_[0]};
	my $count = 0;
	my $sum = 0;
	for my $num(@arr){
		if (defined $num){
			$count++;
			$sum += $num;
		}
	}
	if ($count == 0){
		return 0;
	}
	return $sum/$count;
	
}


my @h3_host_shift = (2,	3,	4,	9,	10,	11,	14,	16,	18,	19,	20,	22,	23,	25,	47,	66,	69,	73,	78,	79,	83,	97,	98,	99,	108,	110,	137,	142,	153,	159,	160,	161,	162,	176,	179,	206,	208,	212,	229,	230,	238,	244,	260,	264,	285,	291,	323,	328,	329,	400,	402,	462,	492,	506,	541);
my @h1_host_shift = (2, 9, 14, 15, 22, 47, 61, 62, 71, 73, 78, 85, 88, 89, 97, 100, 102, 113, 130, 132, 138, 144, 146, 149, 151, 153, 154, 155, 157, 167, 168, 169, 171, 172, 173, 176, 177, 180, 186, 200, 201, 202, 203, 205, 206, 209, 210, 211, 212, 218, 223, 227, 230, 232, 235, 238, 240, 243, 251, 252, 257, 261, 268, 274, 275, 277, 278, 285, 286, 288, 289, 293, 294, 299, 302, 314, 323, 324, 325, 326, 331, 337, 389, 415, 420, 434, 453, 459, 466, 470, 515, 541, 542, 548);
my @n1_host_shift = (3, 5, 8, 12, 13, 14, 16, 20, 26, 29, 34, 40, 41, 42, 43, 46, 47, 51, 52, 53, 59, 64, 66, 67, 69, 70, 71, 72, 74, 75, 76, 78, 79, 80, 81, 82, 83, 85, 93, 95, 99, 101, 105, 111, 114, 116, 136, 149, 157, 189, 195, 200, 206, 210, 211, 214, 220, 221, 222, 223, 232, 241, 250, 257, 258, 263, 264, 267, 273, 274, 285, 287, 288, 289, 309, 311, 329, 339, 340, 341, 351, 354, 355, 365, 367, 369, 382, 386, 388, 390, 393, 394, 396, 427, 430, 432, 434, 451, 454, 455);
my @n2_host_shift = (7, 9, 19, 22, 24, 26, 28, 31, 33, 38, 39, 40, 41, 42, 44, 45, 48, 50, 51, 52, 57, 58, 59, 60, 62, 66, 69, 70, 72, 73, 77, 79, 81, 83, 85, 86, 93, 95, 100, 113, 116, 125, 126, 143, 147, 149, 150, 155, 187, 192, 199, 206, 210, 212, 216, 220, 221, 234, 238, 257, 267, 275, 283, 284, 286, 290, 296, 305, 308, 310, 311, 312, 313, 315, 328, 331, 332, 336, 338, 342, 347, 356, 360, 367, 368, 369, 370, 378, 380, 381, 384, 385, 386, 390, 393, 396, 399, 400, 401, 403, 415, 431, 435, 437, 445, 466);
my @h1_host_shift_001 = (203, 168, 299, 251, 288, 201, 167, 252, 302, 62, 9, 238, 314, 324, 275, 285, 154, 172, 176, 459, 420, 2, 211, 202, 130, 470, 274, 257, 14, 323, 89, 294, 261, 235, 100, 286, 415, 200, 206, 15, 85, 78, 210, 71, 453, 466, 337, 22);
my @h3_host_shift_001 = (16, 108, 229, 244, 79, 73, 83, 161, 260, 20, 9 );
my @n2_host_shift_001 = (386, 384, 381, 328, 83, 70, 81, 192, 51, 147, 125, 283, 41, 286, 77, 72, 378, 331, 126, 155, 50, 62, 338, 369, 60, 315, 216, 399, 396) ;
my @n1_host_shift_001 = (189, 382, 214, 340, 311, 274, 157, 430, 455, 74, 341, 220, 221, 288, 351, 80, 264, 289, 365, 339, 52, 46, 59, 42, 34, 47, 393, 427, 3, 67, 309, 329, 29);

# Caton 1982
my @h1_epitopes = qw(141 142 171 173 175 176 178 179 180 169 172 205 206 209 211 153 156 158 182 186 195 220 237 238 253 286 87 88 90 91 92 132);
# Hensley 2009
my @h1_increased_binding = qw(141 142 171 178 179 180 169 193 205 209 211 156 237 87 88 132 257 175 158 106);
my @n1_epitopes = qw(380 381 382 383 384 385 386 388 389 390 393 397 398 199 200 201 202 223 329 330 332 331 333 336 337 339 340 341 343 344 356 363 364 365 366 367);

# Whiley 
my @h3_epitopes = qw( 138 140 142 147 148 146 149 151 153 154 156 158 159 160 161 162 166 168 184 144 145 171 172 173 174 175 176 179 180 181 202 203 204 205 206 208 209 210 212 213 214 60 61 62 63 64 66 67 69 70 289 291 292 294 295 296 310 313 315 316 320 321 323 324 325 326 327 328 112 118 119 133 137 183 186 187 188 189 190 191 192 193 195 198 217 219 223 224 225 228 229 230 231 232 233 234 235 242 243 244 245 246 254 256 258 260 262 263 264 73 75 78 79 83 91 94 96 97 98 99 102 103 104 107 108 110 125 276 277 278 281 );

# Shih 2007
my @h3_shih_epitopes = qw(66 69 70 137 138 140 142 147 149 151 153 158 159 160 161 162 171 172 173 174 175 176 179 180 188 189 190 202 204 205 206 208 209 212 213 217 223 229 233 242 243 258 260 264 291 292 294 315 323 );
my @n2_epitopes = qw(383 384 385 386 387 389 390 391 392 393 394 396 399 400 401 403 197 198 199 200 221 222 328 329 330 331 332 334 336 338 339 341 342 343 344 346 347 357 358 359  366 367 368 369 370);
my @n1_wan_epitopes = qw(248, 249, 250, 273, 309, 338, 339, 341, 343, 396, 397, 456);
#h1 Huang (antigenic), as is in file Tables_main (Huang + 17, from msa)
my @h1_antigenic = qw( 138 144 145 147 150 158 163 142 170 177 200 203 206 207 208 210 211 52 53 60 288 290 291 294 312 327 111 180 222 226 233 239 241 64 71 86 88 90 97 99 284 );


#h1 Ren, Li, Liu (antigenic) - intersection of two methods
my @h1_antigenic_ren = qw(60 71	88	138	142	144	147	158	204	207	210	222	338);

my @h1_antigenic_Huang_and_host_shift = qw(71	200	203	206	210	211	288	294	);




my @n2_surface = (34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 95, 107, 110, 111, 112, 113, 118, 125, 126, 127, 128, 130, 141, 143, 146, 147, 149, 150, 151, 152, 153, 154, 160, 161, 162, 169, 171, 173, 187, 189, 196, 197, 198, 199, 200, 208, 209, 210, 212, 215, 216, 218, 219, 220, 221, 222, 224, 234, 236, 244, 245, 246, 247, 248, 249, 250, 251, 253, 258, 259, 261, 262, 263, 264, 265, 267, 268, 269, 270, 271, 273, 277, 283, 284, 285, 286, 292, 295, 296, 304, 306, 307, 308, 309, 310, 311, 312, 313, 315, 326, 328, 329, 330, 331, 332, 334, 336, 337, 338, 339, 341, 342, 343, 344, 346, 347, 356, 357, 358, 359, 366, 367, 368, 369, 370, 371, 378, 380, 381, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 396, 399, 400, 401, 402, 403, 413, 414, 415, 416, 417, 430, 431, 432, 433, 434, 435, 437, 450, 451, 452, 453, 455, 456, 457, 459, 461, 463, 464, 465, 466, 468, 469, 470);
my @n1_surface = (34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 86, 88, 89, 90, 93, 94, 95, 111, 118, 126, 127, 128, 136, 141, 143, 146, 147, 148, 149, 150, 151, 152, 154, 162, 163, 165, 172, 174, 189, 191, 199, 200, 201, 202, 209, 210, 211, 214, 215, 217, 218, 220, 221, 222, 223, 224, 226, 236, 237, 247, 248, 249, 250, 251, 252, 255, 258, 260, 261, 263, 264, 265, 266, 267, 268, 269, 271, 273, 274, 275, 279, 285, 286, 287, 288, 290, 296, 297, 298, 304, 306, 307, 308, 309, 311, 312, 313, 314, 326, 328, 329, 330, 331, 332, 333, 334, 335, , , 336, 337, 338, 340, 341, 349, 350, 351, 352, 360, 361, 362, 363, 364, 365, 372, 374, 375, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 389, 391, 392, 393, 394, 395, 398, 406, 407, 408, 413, 414, 415, 416, 417, 426, 427, 429, 430, 431, 432, 433, 434, 435, 437, 439, 450, 451, 452, 454, 455, 456, 457, 461, 463, 464, 465, 467, 468);

my @n1_internal = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,85,87,91,92,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,115,116,117,119,120,121,122,123,124,125,129,130,131,132,133,134,135,137,138,139,140,142,144,145,153,155,156,157,158,159,160,161,164,166,167,168,169,170,172,174,175,176,177,178,179,180,181,182,183,184,185,186,187,189,191,192,193,194,195,196,197,202,203,204,205,206,207,211,212,215,218,224,226,227,228,229,230,231,232,233,234,237,238,239,240,241,242,243,244,245,252,253,255,256,258,261,269,271,275,276,277,279,280,281,282,283,288,290,291,292,293,294,298,299,300,301,302,304,306,310,315,316,317,318,319,320,321,322,323,324,325,327,333,342,345,346,347,348,349,350,351,356,357,358,359,360,361,362,369,370,371,372,373,374,376,379,380,391,393,394,400,401,403,404,405,406,407,408,409,418,419,420,421,422,423,424,425,428,434,439,441,442,443,444,445,446,447,448,449,450,454,459,460,461,463,467,470,471,472);
my @n2_internal = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 87, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 108, 109, 114, 115, 116, 117, 119, 120, 121, 122, 123, 124, 129, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 142, 144, 145, 148, 155, 156, 157, 158, 159, 163, 164, 165, 166, 167, 168, 170, 172, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 188, 190, 191, 192, 193, 194, 195, 201, 202, 203, 204, 205, 206, 207, 211, 213, 214, 217, 223, 225, 226, 227, 228, 229, 230, 231, 232, 233, 235, 237, 238, 239, 240, 241, 242, 243, 252, 254, 255, 256, 257, 260, 266, 272, 274, 275, 276, 278, 279, 280, 281, 282, 287, 288, 289, 290, 291, 293, 294, 297, 298, 299, 300, 301, 302, 303, 305, 314, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 327, 333, 335, 340, 345, 348, 349, 350, 351, 352, 353, 354, 355, 360, 361, 362, 363, 364, 365, 372, 373, 374, 375, 376, 377, 379, 382, 395, 397, 398, 404, 405, 406, 407, 408, 409, 410, 411, 412, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 436, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 454, 458, 460, 462, 467);

#h3 antigentic Steinbruck
my @h3_antigenic  = qw( 138 160 171 223 161 205 233 294 66 153 174 276 140 151 230 278 78 172 212 292 41 91 99 147 202 218 238 241 19 204 69 180 190 209 217 229 246 98 149 159 162 176 213 18 70 188 260 206 242 );

# h3 antigenic change Koel
my @h3_antigenic_koel = (161, 171, 172, 174, 175, 205, 209);
my @h3_antigenic_smith = (138, 160, 153, 161, 149, 159, 162, 140, 147, 171, 204, 180, 205, 209, 174, 172, 176, 213, 175, 206, 212, 294, 66, 70, 223, 190, 217, 229, 233, 246, 188, 260, 292, 98, 276, 78, 91, 99, 69);
my @h3_antigenic_steinbruck  = (138, 160, 171, 223, 161, 205, 233, 294, 66, 153, 174, 276, 140, 151, 278, 230, 78, 172, 212, 292, 41, 91, 99, 147, 202, 218, 238, 241);

#my @h1_pocket_closest = (207,239,166,211,148,111,208,206,234,149,235,196,203,210,151,150,241,204,237,236,240,238,205,209,242,110,197,195,147,202,212,165,152,112,167,233,158,243,244,198,265,159,199,263,201,153,169,262,160);
my @h1_pocket_closest = (207,239,166,111,149,196,203,241,147,168,240,242,110,169,146,197,195,148,206,238,150,204,208,202,167,165,112,243,244,235,198,265,158,199,263,144,205,151,236,159,201,237,152,145,200,143,113,209,245);
my @h3_pocket_closest = (152,241,242,206,238,210,151,114,153,169,165,205,209,42,29,28,211,239,150,207,237,168,113,154,170,115,240,243,212,208,166,43,162,164,161,171,268,175,27,203,204,163,149,155,244,167,172,41,30,148,236,44,156,269);
my @n2_pocket_closest = (151,277,406,150,407,276,278,152,405,226,350,425,291,225,292);
my @n2_bigger_pocket_closest = (292,152,119,276,224,227,406,118,277,371,117,225,407,275,278,372,370,291,223,153,151,228,405,293,226,120,350,425,180,404,181,242,178,133,134,300,241,440,179,441,365,222,150,349);
my @n1_pocket_closest = (402,151,278,152,403,401,277,279,150,227,425,347,292,424,226,228);

my @h1_surface = (28,30,38,39,40,42,52,53,54,56,62,63,64,68,71,83,86,91,92,100,101,103,111,132,136,137,138,141,142,146,148,150,155,156,157,158,171,172,173,176,178,182,184,186,200,201,202,205,206,211,212,221,227,232,235,237,238,252,253,277,278,283,285,287,289,290,292,294,299,303,305,313,324,326,327,339,340,344,350,351,354,358,359,361,362,370,372,373,375,377,381,382,385,386,389,392,396,400,403,404,406,410,412,415,416,425,464,467,470,471,474,476,477,478,484,485,486,488,490,493,497,498,499,501,502,503);
my @h1_internal = (20,22,23,24,25,26,33,35,36,41,43,44,49,50,58,59,66,67,69,72,74,75,76,77,78,79,80,81,82,84,87,93,95,96,97,98,104,105,107,108,109,112,117,118,119,121,122,125,128,133,134,140,143,149,151,152,160,161,163,164,165,166,167,174,177,183,189,190,191,192,193,194,195,196,197,198,199,204,208,213,215,216,217,218,219,226,231,233,241,242,243,245,246,247,248,250,254,256,258,259,260,262,263,264,265,266,267,270,271,272,273,280,282,284,293,295,296,297,298,300,301,302,307,308,309,315,316,319,320,323,325,328,330,331,333,334);

my @h3_surface = (25,37,38,41,43,47,48,49,54,61,62,64,66,71,73,79,94,97,98,99,107,108,110,112,120,120,140,142,144,147,148,149,151,153,156,158,159,160,161,173,174,175,176,178,179,181,183,187,188,189,204,205,206,208,209,214,215,224,228,230,238,240,241,255,256,277,278,279,280,285,287,289,292,294,301,305,307,326,328,329,340,341,352,353,356,361,363,364,372,374,375,376,377,383,384,387,391,394,398,402,403,405,406,414,416,418,427,466,466,480,488,491,492,499,500,501,503,505,506,509,510,513,517,518);
my @h3_internal = (27,29,31,32,33,35,42,44,45,50,52,53,58,59,60,67,68,72,75,77,80,82,83,84,85,86,87,88,89,92,95,100,102,103,104,105,106,113,114,115,118,123,124,125,126,127,128,129,131,132,133,134,136,141,143,146,155,163,164,166,167,168,169,170,177,180,182,186,192,193,194,195,196,197,198,199,200,201,202,207,211,216,218,219,220,221,222,229,231,236,244,245,246,247,248,251,253,257,259,260,261,263,265,266,267,268,269,270,272,273,274,281,282,283,284,286,288,291,297,298,299,302,303,304,310,311,318,319,321,322,325,330,332,333,335,336,338,349,350,351,355,358,359,362,367,368,369,373,385,386,389,393,396,411,420,423,425,426,428,429,430,432,434,435,436,437,438,440,441,444,445,446,448,449,452,453,454,455,456,457,458,459,460,463,464,467,469,471,474,476,477,481,483,485,486,487,489,493,494,497,502,508,511,512,515);


my @h1_leading_kr = (4,11,13,16,52,60,73,74,86,88,90,91,97,99,111,113,128,144,151,156,157,162,169,170,171,172,173,178,182,184,199,202,203,205,206,207,209,210,232,240,261,268,269,283,287,289,290,293,326,361,415,488,489);
my @h1_trailing_kr = (3,6,7,11,52,53,64,89,91,99,101,111,129,137,142,144,148,150,156,157,158,162,165,169,172,178,179,185,186,195,197,199,200,201,203,207,231,236,243,251,253,269,274,278,289,290,293,299,324,331,361,389,390,398,415,422,455,467,470,489,510,514,515,526,562);
my @h3_leading_kr = (27,35,57,82,89,94,107,115,138,153,156,163,165,167,169,172,175,176,177,187,188,190,191,192,195,204,218,221,222,224,225,229,234,249,251,254,258,259,293,294,307,308,310,379,393,407,418,482,484,561);
my @h3_trailing_kr = (19,26,27,29,32,59,65,77,79,80,81,82,83,85,88,89,107,117,120,124,126,128,138,144,160,163,169,170,172,174,182,189,191,192,195,196,203,204,205,206,224,225,226,231,233,234,239,241,246,248,252,254,255,257,258,261,276,280,292,301,305,308,311,323,355,358,393,407,417,418,450,458,482,500,521,532,554,562,579);
my @n2_leading_kr = (18,20,23,30,52,93,143,150,194,197,199,208,216,220,221,249,265,307,308,310,313,328,336,339,344,346,368,369,370,372,381,385,387,390,432);
my @n2_trailing_kr = (2,4,5,9,27,30,40,44,45,50,56,65,77,82,83,120,127,147,148,149,151,155,210,216,220,238,248,251,258,263,265,269,302,307,309,310,312,328,329,334,335,338,339,342,347,372,386,392,400,402,403,414,416,432,433,434,455,464);
my @n1_leading_kr = (15,17,23,34,45,64,70,78,105,173,200,214,222,234,248,249,250,254,270,274,275,287,329,332,336,339,344,352,354,367,369,382,390,396,418,427,430,434,451);
my @n1_trailing_kr = (15,17,21,23,38,39,40,42,45,47,48,52,57,67,68,70,73,77,81,82,83,93,100,101,114,130,147,149,188,200,249,254,259,262,264,267,270,273,275,329,331,340,346,352,364,366,367,390,416,418,419,427,435,452,453,455,462);


my @h3_leading_neva = (18,176,159,175,243,242,162,260,264,161,377,213,202,19,292,69,208,178,531,138,172);
my @h3_trailing_neva = (97,79,236,342,538,13,179,468,289,264,235,173,400);
my @n2_leading_neva = (126,56,249,332,399,264,431,215,43,267,248,220,290,328,313,46,208,432,172,372,69,308,370);
my @n2_trailing_neva = (401,372,155,335,220,339,432,430,44,127,263);

my @h1_leading_neva = (4,113,171,257,178,138,290,52,415,184,223,16,102,86,199,157,162,361,200,111,74,145,221,169,64);
my @h1_trailing_neva = (91,200,232,169,268);
my @n1_leading_neva = (163,263,388,6,149,59,78,14,80,101,427,386,200);
my @n1_trailing_neva = (434,275,15,267,83,287);


my @h3_deps_evolving = qw(10 61 151 161 171 174 245 264 347);
## 90 240 277 179 - only in egg-adapted before 1979
my @h1_jianpeng_evolving = qw(156 169 171 203 206 210 238 90 240 277 179);

my @h1_wenfu_evolving = qw(98 110 157 178 202 203 238 176);


# constructs a hash of observation_vectors from our object
sub get_observation_vectors {
	my $self = shift;
	if ($self->{obs_vectors}){return %{$self->{obs_vectors}};}
	else {
		my %res_hash;
		my %hash = $self ->get_hashes();
		foreach my $ind (keys %hash){
			my @arr = make_observation_vector($hash{$ind}{"x"}, $hash{$ind}{"y"});
			$res_hash{$ind} = \@arr;
		}
		return %res_hash;
	}
}

sub shuffle_observation_vectors {
	my $obs_vectors = $_[0];
	my %shuffled_obs_vectors;
	foreach my $ind (keys %{$obs_vectors}){
		my @arr = shuffle_obsv(\@{$obs_vectors->{$ind}});
		$shuffled_obs_vectors{$ind} = \@arr;
	}
	#print Dumper (\%shuffled_obs_vectors);
	return %shuffled_obs_vectors;
}

# constructs our object given the observation_vectors (shuffled)
sub read_observation_vectors {
	my $self = shift;
	my $obs_vectors = shift;
	my %subs_on_node;
	my %nodes_with_sub;

	foreach my $ind(keys %{$obs_vectors}){
		foreach my $set(@{$obs_vectors->{$ind}}){
		if ($set->[2] == 1){
			my $nodname = $set->[0];
			my $p=Substitution->new();
			$p->position($ind);
			$p->ancestral_allele("ATG");
			$p->derived_allele("ATG");
			$subs_on_node{$nodname}{$ind} = $p;
			if (! exists $nodes_with_sub{$ind}){
					$nodes_with_sub{$ind} = ();
			}
			push (@{$nodes_with_sub{$ind}}, \${$self->{static_hash_of_nodes}{$nodname}}); #вытащить из дерева по имени
			#push (@{$nodes_with_sub{$ind}}, \${$static_hash_of_nodes{$nodname}}); #вытащить из дерева по имени
		#print "TEST1 ".${$static_hash_of_nodes{$nodname}}->get_name()."\n"; # часть имен исчезла, а часть - осталась ОО
		#print "TEST2 ".$nodname."\n";
		}
		}
	}
	return (\%subs_on_node, \%nodes_with_sub);
}




sub myclone {
	my $self = shift;
	my $clone = {
			static_protein => $self->{static_protein},
			static_tree =>  $self->{static_tree},
			static_fasta => $self->{static_fasta},
			static_state  => $self->{static_state},
			static_hash_of_nodes => $self->{static_hash_of_nodes},
			static_distance_hash => $self->{realdata}{"distance_hash"},
			static_background_subs_on_node => $self->{static_background_subs_on_node },
			static_background_nodes_with_sub => $self->{static_background_nodes_with_sub},
			obs_vectors => clone($self->{realdata}{"obs_vectors"}) ; #the only structure which can (and will) be changed
	};
	
	return $clone;
}



sub shuffle_mutator {
	my $self = shift;
	my %obs_vectors = $self->get_observation_vectors();
	my %shuffled_obs_vectors = shuffle_observation_vectors(\%obs_vectors);
	$self->{obs_vectors} = \%shuffled_obs_vectors;
	my @mock_mutmaps = $self->read_observation_vectors(\%shuffled_obs_vectors); 
	$self->{static_subs_on_node} = $mock_mutmaps[0];
	$self->{static_nodes_with_sub} = $mock_mutmaps[1];
	return $self; 
}

# 5.11 for entrenchment_bootstrap_full_selection_vector
sub iterations_gulp {
	my $self = shift;
	my $iterations = shift;
		
	#my $ancestor_nodes = $_[4];
	#my $obs_vector = $_[5];
	#my $norm = $_[6];
		
	my $realdata = $self->{realdata};
	my $maxbin = $realdata->{"maxbin"};
	my $ancestor_nodes = $realdata->{"ancestor_nodes"};
	#my $obs_vectors = $realdata->{"obs_vectors"};

	my $mock_mutmap = myclone(); # 25.07 create a new object for shuffling
	my @simulated_hists;
	
	for (my $i = 1; $i <= $iterations; $i++){
		$mock_mutmap->shuffle_mutator(); #this method shuffles observation vectors and sets new $static_nodes.. and static_subs..
		my %hash;
		# >new iteration string and all the corresponding data  are printed inside this sub:
		my %prehash = $mock_mutmap->depth_groups_entrenchment_optimized_selection_alldepths(1,0,$ancestor_nodes, "overwrite"); #bin, restriction (NOT USED), ancestor_nodes, should I overwrite static hash?

		foreach my $bin(1..$maxbin){
				foreach my $site_node(keys %prehash){
					$hash{$bin}[1] += $prehash{$site_node}{$bin}[1];
					$hash{$bin}[0] += $prehash{$site_node}{$bin}[0];				
			}
		}
		# 21.12 we do not need any norm here since we normalize values in concat_and_divide - separately for different maxdepths
		push @simulated_hists, \%hash;
		#%static_ring_hash = ();
		#%static_subtree_info = ();
	}

	store \@simulated_hists, $self->{static_output_base}.$self->{static_protein}."_gulpselector_vector_alldepths_stored";
	
	my $arref = retrieve($self->{static_output_base}.$self->{static_protein}."_gulpselector_vector_alldepths_stored");
	my $csvfile = $self->{static_output_base}.$self->{static_protein}."_gulpselector_vector_alldepths.csv";
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
	
}

# 28.12 compute norm for given restriction and/or group
sub compute_norm {
	my $self = shift;
	my $restriction = $_[0];
	my @group;
	if ($_[1]){
		@group = @{$_[1]};
	}
	else {
		@group = (1..565);
	}
	my %group_hash;
	foreach my $ind(@group){
		$group_hash{$ind} = 1;
	}
	
#	$real_data = lock_retrieve("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_realdata") or die "Cannot retrieve real_data";
	my $realdata = $self->{realdata};
	@obshash_restriction = map { /^obs_hash(.*)/ ? $1 : () } (keys $realdata);
	unless(defined $obshash_restriction[0] && $obshash_restriction[0] <= $restriction ){
			die "realdata restriction is bigger than compute_norm restriction: ".$obshash_restriction[0]." > $restriction \n";
	}
	my $obs_hash = $realdata->{"obs_hash".$obshash_restriction[0]};
	my $subtree_info = $realdata->{"subtree_info"};
	
	my $norm;
	
	foreach my $site_node(keys %{$obs_hash}){
		my ($site, $node_name) = split(/_/, $site_node);
		if ($group_hash{$site}){
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			foreach my $bin(keys %{$obs_hash->{$site_node}}){
			
				if ($maxdepth > $restriction){
					$norm += $obs_hash->{$site_node}->{$bin}->[0];
				}

			}
		}
	}
	return $norm;
}

#26.02 N groups

sub print_nodes_in_analysis {
	my $self = shift;
	my $prot = $self->{static_protein};
	my $restriction = $_[0];
	my @groups = @{$_[1]};
	my @names = @{$_[2]};
	my $subtract_tallest = $self->{static_subtract_tallest};
	$self->set_distance_matrix();
	#my %matrix = incidence_matrix(); #!
	my $i = 0;
	foreach my $group(@groups){
		print $names[$i]."\n";
		my %obs_hash = $self->nodeselector(1,$restriction,$subtract_tallest, $group, $names[$i]); #bin, restriction, subtract-tallest

		%static_ring_hash = ();
		%static_subtree_info = ();
		$i++;
	}
}



#5.11 for entrenchment_bootstrap_full_selection_vector
# analyze real data, prepare bkg_mutmap, ancestor_nodes and obs_vector
#28.12 hash is pruned: we do not keep mutation info if its maxdepth<50
#1.08 hash is not necessarily pruned - you can set restriction to 0 to get complete data
sub prepare_real_data {
	my $self = shift;
	my $restriction = shift;
	unless(defined $restriction) { $restriction = 50; }
	$self -> set_distance_matrix();
	my %matrix = $self->incidence_matrix(); 
	$self -> print_incidence_matrix(\%matrix);
	# used depth_groups_entrenchment_optimized_selector_alldepths but changed it for depth_groups_entrenchment_optimized_selector_alldepths_2, because the latter
	# keeps in obs_hash info about site index as well as about node name
	my %full_obs_hash = $self -> depth_groups_entrenchment_optimized_selector_alldepths_2(1, $restriction); # bin size
	my %ancestor_nodes;
	foreach my $ancnode(keys %full_obs_hash){
	my @splitter = split(/_/, $ancnode);
		$ancestor_nodes{$splitter[-1]} = 1;
	}	
	
	my $restricted_norm;
	my %restricted_obs_hash;
	my $maxbin = 0;
	
	foreach my $site_node(keys %full_obs_hash){
		my ($site, $node_name) = split(/_/, $site_node);
		my $maxdepth = $self -> {static_subtree_info}{$node_name}{$site}{"maxdepth"};
		foreach my $bin(keys %{$full_obs_hash{$site_node}}){
			$maxbin = max($bin, $maxbin);
			if ($maxdepth > $restriction){
				$restricted_norm += $full_obs_hash{$site_node}{$bin}[0];
				$restricted_obs_hash{$site_node}{$bin}[0] = $full_obs_hash{$site_node}{$bin}[0];
				$restricted_obs_hash{$site_node}{$bin}[1] = $full_obs_hash{$site_node}{$bin}[1];
			}
		}
#		print "MAXBIN $maxbin\n";
	}
	
#	print " NORM50 $norm50  NORM100 $norm100  NORM150 $norm150 \n";
	my %obs_vectors = $self ->get_observation_vectors();

	my %realdata = (
		norm."restriction" => $restricted_norm,
		restriction = $restriction,
		maxbin => $maxbin,
		ancestor_nodes => \%ancestor_nodes,
		obs_vectors => \%obs_vectors,
		bkg_subs_on_node => $self -> {static_background_subs_on_node},
		bkg_nodes_with_sub => $self -> {static_background_nodes_with_sub},
		distance_hash => $self -> {static_distance_hash},
		hash_of_nodes => $self -> {static_hash_of_nodes},
		subtree_info => $self -> {static_subtree_info},
		obs_hash."restriction" => \%restricted_obs_hash,
	);
	
	my $realdatapath = $self->{static_output_base};
	$realdatapath = File::Spec->catfile($realdatapath, $prot."_".state_tag($self->{state})."_realdata");
	print "Saving real_data to $realdatapath\n";

	store \%realdata, $realdatapath;
	#$self->{static_ring_hash} = ();
	#$self->{static_subtree_info} = ();
	
	#%static_subs_on_node = ();
	#%static_background_subs_on_node = ();
	#%static_nodes_with_sub = ();
	#%static_background_nodes_with_sub = ();
}


	
	
sub select_ancestor_nodes {
	my $self = shift;
	my $restriction = $_[0];
	my @group = @{$_[1]};
	my %group_hash;
	foreach my $site(@group){
		$group_hash{$site} = 1;
	}
	my $realdata = $self->{realdata};
	my $maxbin = $realdata->{"maxbin"};
	my $ancestor_nodes = $realdata->{"ancestor_nodes"};
	@obshash_restriction = map { /^obs_hash(.*)/ ? $1 : () } (keys $realdata);
	unless(defined $obshash_restriction[0] && $obshash_restriction[0] <= $restriction ){
			die "realdata restriction is bigger than select_ancestor_nodes restriction: ".$obshash_restriction[0]." > $restriction \n";
	}
	my $obs_hash = $realdata->{"obs_hash".$obshash_restriction[0]};
	my $subtree_info = $realdata->{"subtree_info"};
	my %group_nodes;
	foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = split(/_/, $site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
				if ($maxdepth > $restriction && $group_hash{$site}){
					$group_nodes{$node_name} = 1;
					#print "group_node ".$node_name."\n";
				}
		}
		my $count = scalar keys %group_nodes;
	#print "Total $count\n";
	return %group_nodes;
}	
	
	

	
	
sub prepare_groups {
	my @groups = @{$_[0]};
	my @final_groups;
	foreach my $pregroup (@groups){
		my %ghash;
		my @group;
		foreach my $s(@{$pregroup}){
			$ghash{$s} = 1;
			push @group, $s;
		}
		my @complement;
		foreach my $s(1..565){
			if (!$ghash{$s}){
				push @complement, $s;
			}
		}
		push @final_groups, \@group;
		push @final_groups, \@complement;
	}
	# both group and complement must contain only variable sites:
	#my %variable_sites;
	#my @final_groups;
	#my %nodes_with_sub = $realdata->{"nodes_with_sub"};
	#foreach my $s(1..565){
	#	if ($nodes_with_sub{$s}){
	#		$variable_sites{$s} = 1;
	#	}
	#}
	#foreach my $pregroup (@groups){
	#	my %ghash;
	#	my @group;
	#	foreach my $s(@{$pregroup}){
	#		if ($variable_sites{$s}){
	#			$ghash{$s} = 1;
	#			push @group, $s;
	#		}
	#	}
	#	my @complement;
	#	foreach my $s(keys %variable_sites){
	#		if (!$ghash{$s}){
	#			push @complement, $s;
	#		}
	#	}
	#	push @final_groups, \@group;
	#	push @final_groups, \@complement;
	#}
	return @final_groups;
}	

sub prepare_groups_and_names {
	my @pregroups = @{$_[0]};
	my @pregroup_names = @{$_[1]};
	
	my @group_names;
	
	foreach my $group_name(@pregroup_names){
		push @group_names, $group_name;
		push @group_names, $group_name."_complement";
	}
	
	my @groups = prepare_groups(\@pregroups);

	#my @all_variable_sites;
	#my %nodes_with_sub = $realdata->{"nodes_with_sub"};
	#foreach my $s(1..565){
	#	if ($nodes_with_sub{$s}){
	#		print " pushed $s\n";
	#		push @all_variable_sites, $s;
	#	}
	#}
	
	my @all_sites = (1..565);
	
	push @groups, \@all_sites;
	push @group_names, "all";
	
	return (\@groups, \@group_names);
}

#16.02 complement to leading is trailing, and vice versa
sub prepare_lt_groups_and_names {
	my $prot = $_[0];
	my @group_names;
	push @group_names, "pure_leading";
	push @group_names, "trailing_as_complement";
	push @group_names, "pure_trailing";
	push @group_names, "leading_as_complement";	

	my @groups = prepare_lt_groups($prot);
	
	my @all_sites = (1..565);
	
	push @groups, \@all_sites;
	push @group_names, "all";
	
	return (\@groups, \@group_names);
}

#16.02 complement to leading is trailing, and vice versa
sub prepare_lt_groups {		
		my $prot = $_[0];
		my @leading;
		my @trailing;
		my @final_groups;
		switch ($prot) {
			case "h1" {	@leading = @h1_leading_kr;
						@trailing = @h1_trailing_kr;
					  }
			case "h3" {	@leading = @h3_leading_kr;
						@trailing = @h3_trailing_kr;
					  }
			case "n1" {	@leading = @n1_leading_kr;
						@trailing = @n1_trailing_kr;
					  }
			case "n2" {	@leading = @n2_leading_kr;
						@trailing = @n2_trailing_kr;
					  }	
			else {	print "Panic in prepare_lt_groups! No such protein $prot!\n";
					die;
				 }			
		}
		my %trailing = map{$_ =>1} @trailing;
		my %leading = map{$_ =>1} @leading;
		
		my @pure_leading = grep(!defined $trailing{$_}, @leading);
		my @pure_trailing = grep(!defined $leading{$_}, @trailing);
	
		push @final_groups, \@pure_leading;
		push @final_groups, \@pure_trailing;
		push @final_groups, \@pure_trailing;
		push @final_groups, \@pure_leading;
		
		return @final_groups;
}




#27.01 writes to files immediately, does not hold hash of iterations in memory	
	
sub concat_and_divide_simult {
	my $self = shift;
	my $prot = $self->{static_protein};
	my @maxdepths = @{$_[0]};
	my @groups = @{$_[1]};
	my @group_names = @{$_[2]};
	my $subtract_maxpath = $self->{static_subtract_tallest};
	my $dir = $self->{static_output_base};
	my $nodecount_file = File::Spec->catfile($dir, $prot."_nodecount");
	open NODECOUNT, ">$nodecount_file" or die "Cannot create $nodecount_file";
	
	my $realdata = $self->{realdata};
	my %hash;
	my $iteration_number = 1;
	my $maxbin = $realdata->{"maxbin"};
	
	my %norms;
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			$norms{$md}[$group_number] = $self->compute_norm($md, $groups[$group_number]);
		}
		
	}
	
	my %group_hashes;
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			my %node_hash = $self->select_ancestor_nodes($md, \@{$groups[$group_number]});
			foreach my $node_name (keys %node_hash){
				$group_hashes{$md}[$group_number]{$node_name} = $node_hash{$node_name};
			}
		}
		
	}
	
	
	
	#foreach my $gulp(1..$gulps){
	#	my $gulp_filename = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_for_enrichment_".$tag.$gulp;
	my $dirname = File::Spec->catdir($dir, $prot);
	
	#if ($syn){
	#	#$dirname = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."_syn/";
	#	 $dirname = File::Spec->catdir($dirname, $prot."_syn");
	#}
	#else {
	#	#$dirname = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/".$prot."/";
	#	$dirname = File::Spec->catdir($dirname, $prot);
	#}
	opendir(DH, $dirname);
	my @files = readdir(DH);
	closedir(DH);
	
	my %filehandles;

	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			local *FILE;
			my $csvfile =  File::Spec->catfile($dir, $prot."_gulpselector_vector_".$md."_".$group_names[$group_number]."_".$tag.".csv");
			open FILE, ">$csvfile" or die "Cannot create $csvfile";
			
			FILE->autoflush(1);
			$filehandles{$md}{$group_number} = *FILE;
		}
		
	}

	
	
	
	foreach my $gulp_filename(@files){
	next if (-d $gulp_filename);
	my $fullpath = File::Spec -> catfile($dirname, $gulp_filename);
	open GULP, "<$fullpath" or die "Cannot open $fullpath";
	
	my $site; # careful	
	my $node_name;
		while (<GULP>){

			my $max_depth;
		#5.02	if ($_ =~ /^>/){
		#5.02		#$iteration_number++;
		#5.02		my $str = <GULP>;
		#5.02		#print "str ".$str."\n";
		#5.02		my @str_array = split(/\s+/, $str);
		#5.02		my $site = $str_array[1];
		#5.02		my $node_name = $str_array[3];
		#5.02		$max_depth = $str_array[5];			
		#5.02		#print $site." site,".$node_name." node,".$max_depth." maxdepth\n";
		#5.02	}
			my @str_array;
			my $str = <GULP>;
			
			my %sums;
			my %hash;
			
			my $test_obs_summ;
			my $test_exp_summ;
			
			while ($str =~ /^[^>]/){ 

			
				if ($str =~ /^site/){
				
			if ($test_obs_summ - $test_exp_summ < 0.0001 && -$test_obs_summ + $test_exp_summ < 0.0001){
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
					
					## 15.02 iterations: count number of nodes in each group
					foreach my $md(@maxdepths){
						if ($max_depth > $md){
							foreach my $group_number(0..scalar @groups-1){
								if ($group_hashes{$md}[$group_number]{$node_name}){
									print NODECOUNT "maxdepth $md group ".$group_names[$group_number]." node $node_name\n";
								}
							}
						}	
					}
					##
					
					
					
					#print $site." site,".$node_name." node,".$max_depth." maxdepth\n";
					#5.02 my $str = <GULP>; 
					$str = <GULP>; #5.02
				}
				@str_array = split(/,/, $str);
				
								$test_obs_summ += $str_array[1];
								$test_exp_summ += $str_array[2];
				
								#print " Maxdepth $max_depth itnum $iteration_number bin ".$str_array[0]." exp ".$str_array[2]." obs ".$str_array[1]." \n";
								foreach my $md(@maxdepths){
									if ($max_depth > $md){
										foreach my $group_number(0..scalar @groups-1){
											if ($group_hashes{$md}[$group_number]{$node_name}){
											#print "group number $group_number md $md node name $node_name\n";
												$sums{$md}[$group_number] += $str_array[1];
												$hash{$md}[$group_number]{$str_array[0]}[1] += $str_array[2];
												$hash{$md}[$group_number]{$str_array[0]}[0] += $str_array[1];
												#print $hash{$md}[$group_number]{$str_array[0]}[0]." obs integral\n";
											}
										}
									}
								}

				
				$str = <GULP>;
				
			}
			

			# maxbins are different for every iteration. Find maximum and use it.

			$maxbin = max($maxbin, $str_array[0]);
#print "maxbin $maxbin, iteration number $iteration_number\n";	
#print "sum50 $sum50 sum100 $sum100 sum150 $sum150 norm 50 $norm50 norm 100 $norm100 norm 150 $norm150\n";		
			
			foreach my $md(@maxdepths){ 
				foreach my $group_number(0..scalar @groups-1){
					#print "maxdepth $md group number $group_number \n";
					if ($sums{$md}[$group_number] == 0){
						foreach my $bin(1..$maxbin){
							$hash{$md}[$group_number]{$bin}[0] = "NA";
							$hash{$md}[$group_number]{$bin}[1] = "NA";
						}
					}
					else {
						foreach my $bin(1..$maxbin){
							#print "in hash: ".$hash{$md}[$group_number]{$bin}[0]."\n";
							#print "norm ".$norms{$md}[$group_number]."\n";
							#print "sum ".$sums{$md}[$group_number]."\n";
							$hash{$md}[$group_number]{$bin}[0] = $hash{$md}[$group_number]{$bin}[0]*$norms{$md}[$group_number]/$sums{$md}[$group_number];
							$hash{$md}[$group_number]{$bin}[1] = $hash{$md}[$group_number]{$bin}[1]*$norms{$md}[$group_number]/$sums{$md}[$group_number];
						}
					}
					
					my $filehandle = $filehandles{$md}{$group_number};
					print "going to print something\n";
					foreach my $bin(1..$maxbin){
						print $filehandle $hash{$md}[$group_number]{$bin}[0].",".$hash{$md}[$group_number]{$bin}[1].",";
					}
					print $filehandle "\n";
				
				
				}
			}
			

			
			$iteration_number++;
			print $iteration_number."\n";
		}
		
		close GULP;
		
	}
	
	close NODECOUNT;
	
	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
					my $filehandle = $filehandles{$md}{$group_number};
					close $filehandle;
		}
		
	}
	
	
	
	
#	# write to file: out of this cycle
#			open CSV50, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_gulpselector_vector_50_".$tag.".csv";
#			open CSV100, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_gulpselector_vector_100_".$tag.".csv";
#			open CSV150, ">/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_gulpselector_vector_150_".$tag.".csv";
#			# normalization: in cycle or also out of it, if sums are kept in separate hash
#			foreach my $i (1..$iteration_number-1){
#				foreach my $bin(1..$maxbin){
#					print CSV50 $hash{50}[$i]{$bin}[0].",".$hash{50}[$i]{$bin}[1].",";
#					print CSV100 $hash{100}[$i]{$bin}[0].",".$hash{100}[$i]{$bin}[1].",";
#					print CSV150 $hash{150}[$i]{$bin}[0].",".$hash{150}[$i]{$bin}[1].",";
#				}
#				print CSV50 "\n";
#				print CSV100 "\n";
#				print CSV150 "\n";
#			}	
#			close CSV50;
#			close CSV100;
#			close CSV150;
}	
	
	
# 19.12 after concat_and_divide	
sub count_pvalues{	


my $prot = $_[0];
my $tag = $_[1];
my @restriction_levels = @{$_[2]};
my @groups = @{$_[3]};
my @group_names = @{$_[4]};
my $dir = $_[5];


#my $countfile = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/";
#if ($locally){
#	$countfile = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/";
#}

my $countfile = File::Spec->catfile($dir, $prot."_count_".$tag);
open COUNTER, ">$countfile" or die "Cannot create $countfile";
COUNTER->autoflush(1);

my $realdata = get_real_data();

my $maxbin = $realdata->{"maxbin"}; 
#print "before cycle\n";
my $obs_hash = $realdata->{"obs_hash50"};
my $subtree_info = $realdata->{"subtree_info"};
for my $restriction(@restriction_levels){
print "level $restriction\n";

	# only for @all@
	 my $group_number = scalar @groups - 1;
	 print " groups ".scalar @groups - 1;
	my %group_hash;
	print " size ".scalar @{$groups[$group_number]};
	foreach my $site(@{$groups[$group_number]}){
		print "test-1 $site\n";
		$group_hash{$site} = 1;
	}
	my %obs_hash_restricted;
	my $norm_restricted;
## copypaste from prepare: create restricted hash	
# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
	foreach my $site_node(keys %{$obs_hash}){
		my ($site, $node_name) = split(/_/, $site_node);
		my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
		if ($maxdepth > $restriction && $group_hash{$site}){ #15.02 moved here
		foreach my $bin(keys %{$obs_hash->{$site_node}}){
			$maxbin = max($bin, $maxbin);
			#if ($maxdepth > $restriction && $group_hash{$site}){ #15.02 moved from here
				$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
				$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
				$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
			}
	
		}
#		print "MAXBIN $maxbin\n";
	}
	my $count = scalar keys %obs_hash_restricted;
	print COUNTER "Total for $restriction all $count\n";

	
	
	
## end of copypaste	
	
	#my $file = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/";
	#if ($locally){
	#	$file = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/";
	#}
	#$file = $file.$prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]."_".$tag;
	
	my $file = File::Spec->catfile($dir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]."_".$tag);
	open FILE, ">$file" or die "Cannot create $file";
			
	my %histhash;
	foreach my $site_node(keys %obs_hash_restricted){
		foreach my $bin (keys %{$obs_hash_restricted{$site_node}}){
			$histhash{$bin}[0] += $obs_hash_restricted{$site_node}{$bin}[0];
			$histhash{$bin}[1] += $obs_hash_restricted{$site_node}{$bin}[1];
		}
	}
	print FILE "bin\tobs\texp\n";
	my @sorted_bins = sort { $a <=> $b } keys %histhash;
	foreach my $bin (@sorted_bins){
		print FILE $bin."\t".$histhash{$bin}[0]."\t".$histhash{$bin}[1]."\n";
	}
		
		
		my %flat_obs_hash;
		my %flat_exp_hash;
#print " going to flat hash\n";
		foreach my $bin(1..$maxbin){
			foreach my $node (keys %obs_hash_restricted){
				$flat_obs_hash{$bin} += $obs_hash_restricted{$node}{$bin}[0]; 
				$flat_exp_hash{$bin} += $obs_hash_restricted{$node}{$bin}[1]; 
			}
		}
#print " computing hist emdian \n"	;
		my $obs_median = hist_median_for_hash(\%flat_obs_hash);
		my $exp_median = hist_median_for_hash(\%flat_exp_hash);
		my $obs_mean = hist_mean_for_hash(\%flat_obs_hash); # 18.03 - added the same statistics based on histogram mean (instead of median)
		my $exp_mean = hist_mean_for_hash(\%flat_exp_hash);
		
		print FILE "\n observed median: $obs_median\n";
		print FILE "\n poisson expected median: $exp_median\n";
		print FILE "\n observed mean: $obs_mean\n";
		print FILE "\n poisson expected mean: $exp_mean\n";
		
		my $pval_epi;
		my $pval_env;
		my $pval_epi_for_mean;
		my $pval_env_for_mean;

	
#print "going to read input file\n";	
		#my $csvfile = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/";
		#if ($locally){
		#	$csvfile = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/";
		#}
		#$csvfile = $csvfile.$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number]."_".$tag.".csv";
		
		if ($obs_mean ne "NaN"){
		
		my $csvfile = File::Spec->catfile($dir, $prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number]."_".$tag.".csv");
		open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
		my $iteration = 0;
		my @hist_obs;
		my @hist_exp;
		my @array_obs_minus_exp;
		while(<CSVFILE>){
			my %boot_obs_hash;
			my %boot_exp_hash;
			my @splitter = split(/,/, $_);
			my @diff_10bin_array;
			for (my $i = 0; $i < scalar @splitter; $i++){
				my $bin = ($i/2)+1;
				my $obs = $splitter[$i];
				$boot_obs_hash{$bin} = $splitter[$i];
				#print " i $i bin $bin value  $splitter[$i]\n";
				$hist_obs[$bin] += $splitter[$i];
				
				$i++;
				my $exp = $splitter[$i];
				$boot_exp_hash{$bin} = $splitter[$i];
				$hist_exp[$bin] += $splitter[$i];
				$diff_10bin_array[int($bin/10)] += $obs-$exp;
				#push @{$array_obs_minus_exp[$bin]}, $obs-$exp;
			}
			
			for (my $bin10 = 0; $bin10 < scalar @diff_10bin_array; $bin10++){
				push @{$array_obs_minus_exp[$bin10]}, $diff_10bin_array[$bin10];
			}
			
			
			
			
			my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash);
			my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash);
			my $boot_obs_mean = hist_mean_for_hash(\%boot_obs_hash);
			my $boot_exp_mean = hist_mean_for_hash(\%boot_exp_hash);
			print FILE "\n boot obs median: $boot_obs_median boot exp median: $boot_exp_median \n";
			print FILE "\n boot obs mean: $boot_obs_mean boot exp mean: $boot_exp_mean \n";
			if ($boot_obs_median - $boot_exp_median >= $obs_median - $exp_median){
				$pval_env += 1;
			}
			if ($boot_obs_median - $boot_exp_median <= $obs_median - $exp_median){
				$pval_epi += 1;
			}
			if ($boot_obs_mean - $boot_exp_mean >= $obs_mean - $exp_mean){
				$pval_env_for_mean += 1;
			}
			if ($boot_obs_mean - $boot_exp_mean <= $obs_mean - $exp_mean){
				$pval_epi_for_mean += 1;
			}
			$iteration++;
		}
		

			for (my $j = 0; $j < scalar @hist_obs; $j++){
				my $mean_obs = $hist_obs[$j]/$iteration;
				my $mean_exp = $hist_exp[$j]/$iteration;
				my $stat_obs = Statistics::Descriptive::Full->new();
				$stat_obs->add_data(\@{$array_obs_minus_exp[$j]});
				print FILE "bin $j mean_boot_obs $mean_obs mean_boot_exp $mean_exp diff_percentile_5 ".$stat_obs->percentile(5)." diff_percentile_95 ".$stat_obs->percentile(95).".\n";
			}

		
		
		close CSVFILE;
		
		print FILE "- pvalue_epistasis  pvalue_environment\n";
		print FILE "median_stat ".($pval_epi/$iteration)." ".($pval_env/$iteration)."\n";
		print FILE "mean_stat ".($pval_epi_for_mean/$iteration)." ".($pval_env_for_mean/$iteration)."\n";

		close FILE;	
		
		}
		else {
			print FILE "hist sum is 0";
			close FILE;
		}
	

	
	
	## now for the groups (all but "all")
	for (my $group_number = 0; $group_number < scalar @groups - 1; $group_number++){ # only for even

	## group
	my %group_hash;
	
	foreach my $site(@{$groups[$group_number]}){
		$group_hash{$site} = 1;
	}
	my %obs_hash_restricted;
	my $norm_restricted;
## copypaste from prepare: create restricted hash	
# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
	foreach my $site_node(keys %{$obs_hash}){
		my ($site, $node_name) = split(/_/, $site_node);
		my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
		foreach my $bin(keys %{$obs_hash->{$site_node}}){
			$maxbin = max($bin, $maxbin);
			if ($maxdepth > $restriction && $group_hash{$site}){
				$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
				$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
				$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
			}
		}
#		print "MAXBIN $maxbin\n";
	}
	
## end of copypaste	
	my $count = scalar keys %obs_hash_restricted;
	print  COUNTER "$restriction ".$group_names[$group_number]." group $count "; 
		
		
		#my $file = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/";
		#if ($locally){
		#	$file = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/";
		#}
		#$file = $csvfile.$prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]."_".$tag;
		
		my $file = File::Spec->catfile($dir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]."_".$tag);
		open FILE, ">$file";
		
		
		#copypaste from all
		my %histhash;
		foreach my $site_node(keys %obs_hash_restricted){
			foreach my $bin (keys %{$obs_hash_restricted{$site_node}}){
				$histhash{$bin}[0] += $obs_hash_restricted{$site_node}{$bin}[0];
				$histhash{$bin}[1] += $obs_hash_restricted{$site_node}{$bin}[1];
			}
		}
		print FILE "bin\tobs\texp\n";
		my @sorted_bins = sort { $a <=> $b } keys %histhash;
		foreach my $bin (@sorted_bins){
			print FILE $bin."\t".$histhash{$bin}[0]."\t".$histhash{$bin}[1]."\n";
		}
		## end of copypaste
		
		
		my %flat_obs_hash;
		my %flat_exp_hash;
#print " going to flat hash\n";
		foreach my $bin(1..$maxbin){
		foreach my $node (keys %obs_hash_restricted){
			$flat_obs_hash{$bin} += $obs_hash_restricted{$node}{$bin}[0];
			$flat_exp_hash{$bin} += $obs_hash_restricted{$node}{$bin}[1];
		}
		}
#print " computing hist emdian \n"	;
		my $obs_median = hist_median_for_hash(\%flat_obs_hash);
		my $exp_median = hist_median_for_hash(\%flat_exp_hash);
		my $obs_mean = hist_mean_for_hash(\%flat_obs_hash);
		my $exp_mean = hist_mean_for_hash(\%flat_exp_hash);
		
		print FILE "\n observed median: $obs_median expected poisson median $exp_median observed mean: $obs_mean expected poisson mean $exp_mean\n";
	 if($obs_mean eq "NaN" || $exp_mean eq "NaN") {
		print FILE " hist sum is 0";
		close FILE;
		$group_number++; # to skip complement for this group
		next;
	 }
#print "going to read input file\n";		
		
		#my $csvfile = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/";
		#if ($locally){
		#	$csvfile = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/";
		#}
		#$csvfile = $csvfile.$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number]."_".$tag.".csv";
		
		my $csvfile = File::Spec->catfile($dir,$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number]."_".$tag.".csv");
		open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
		my $iteration = 0;
		my @group_boot_medians;
		my @group_boot_means;
		my @hist_obs;
		my @hist_exp;
		my @array_gbo_minus_gbe;
		while(<CSVFILE>){
			my %boot_obs_hash;
			my %boot_exp_hash;
			my @splitter = split(/,/, $_);
			
			## copypaste 5.02
			my @gbo_minus_gbe_10bin_array;
			for (my $i = 0; $i < scalar @splitter; $i++){
				my $bin = ($i/2)+1;
				my $obs = $splitter[$i];
				$boot_obs_hash{$bin} = $splitter[$i];
				#print " i $i bin $bin value  $splitter[$i]\n";
				$hist_obs[$bin] += $splitter[$i];
				
				$i++;
				my $exp = $splitter[$i];
				$boot_exp_hash{$bin} = $splitter[$i];
				$hist_exp[$bin] += $splitter[$i];
				$gbo_minus_gbe_10bin_array[int($bin/10)] += $obs-$exp;
				#push @{$array_obs_minus_exp[$bin]}, $obs-$exp;
			}
			
			for (my $bin10 = 0; $bin10 < scalar @gbo_minus_gbe_10bin_array; $bin10++){
				push @{$array_gbo_minus_gbe[$bin10]}, $gbo_minus_gbe_10bin_array[$bin10];
			}
			## end of copypaste
			
			
			## 5.02 redundant
			#for (my $i = 0; $i < scalar @splitter; $i++){
			#	my $bin = ($i/2)+1;
			#	$boot_obs_hash{$bin} = $splitter[$i];
			#	$i++;
			#	$boot_exp_hash{$bin} = $splitter[$i];
			#}
			
			my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash);
			my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash);
			my $boot_obs_mean = hist_mean_for_hash(\%boot_obs_hash);
			my $boot_exp_mean = hist_mean_for_hash(\%boot_exp_hash);
			
			$group_boot_medians[$iteration][0] = $boot_obs_median;
			$group_boot_medians[$iteration][1] = $boot_exp_median;
			$group_boot_means[$iteration][0] = $boot_obs_mean;
			$group_boot_means[$iteration][1] = $boot_exp_mean;			
			print FILE "\n boot obs median: $boot_obs_median boot obs mean: $boot_obs_mean\n";
			
			$iteration++;
		}
		close CSVFILE;	
		
##complement 

my %complement_hash;
$group_number++;
my $maxbin = 0;
	foreach my $site(@{$groups[$group_number]}){ #next
		$complement_hash{$site} = 1;
	}
	my %complement_obs_hash_restricted;
	my $complement_norm_restricted;
## copypaste from prepare: create restricted hash	
# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
	foreach my $site_node(keys %{$obs_hash}){
		my ($site, $node_name) = split(/_/, $site_node);
		my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
		foreach my $bin(keys %{$obs_hash->{$site_node}}){
			$maxbin = max($bin, $maxbin);
			if ($maxdepth > $restriction && $complement_hash{$site}){
				$complement_norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
				$complement_obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
				$complement_obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
			}
		}
#		print "MAXBIN $maxbin\n";
	}
	
## end of copypaste	
	
		my $count = scalar keys %complement_obs_hash_restricted;
		print  COUNTER " $restriction ".$group_names[$group_number]." complement $count\n"; 
	
		#my $file = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/";
		#if ($locally){
		#	$file = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/";
		#}
		#$file = $csvfile.$prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]."_".$tag;
		
		my $file = File::Spec->catfile($dir,$prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]."_".$tag);
		open FILE, ">$file" or die "Cannot create $file";
		my %complement_flat_obs_hash;
		my %complement_flat_exp_hash;
#print " going to flat hash\n";
		foreach my $bin(1..$maxbin){
			foreach my $node (keys %complement_obs_hash_restricted){
				$complement_flat_obs_hash{$bin} += $complement_obs_hash_restricted{$node}{$bin}[0];
				$complement_flat_exp_hash{$bin} += $complement_obs_hash_restricted{$node}{$bin}[1];
			}
		}
#print " computing hist emdian \n"	;
		my $complement_obs_median = hist_median_for_hash(\%complement_flat_obs_hash);
		my $complement_exp_median = hist_median_for_hash(\%complement_flat_exp_hash);
		my $complement_obs_mean = hist_mean_for_hash(\%complement_flat_obs_hash);
		my $complement_exp_mean = hist_mean_for_hash(\%complement_flat_exp_hash);
		
		print FILE "\n observed median: $complement_obs_median expected median: $complement_exp_median observed mean: $complement_obs_mean expected mean: $complement_exp_mean\n";
		if($complement_obs_mean eq "NaN" || $complement_exp_mean eq "NaN") {
			print FILE " hist sum is 0";
			close FILE;
			next;
		}
#print "going to read input file\n";		
		#my $csvfile = "/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/";
		#if ($locally){
		#	$csvfile = "/cygdrive/c/Users/weidewind/Documents/CMD/Coevolution/Influenza/perlOutput/epi_or_env_december_2015/";
		#}
		#$csvfile = $csvfile.$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number]."_".$tag.".csv";
		
		my $csvfile = File::Spec->catfile($dir,$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number]."_".$tag.".csv");
		open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
		my $iteration = 0;
		my @complement_boot_medians;
		my @complement_boot_means;
		my @hist_compl_obs;
		my @hist_compl_exp;
		my @array_cbo_minus_cbe;
		while(<CSVFILE>){
			my %boot_obs_hash;
			my %boot_exp_hash;
			my @splitter = split(/,/, $_);
			
			## copypaste 5.02
			my @cbo_minus_cbe_10bin_array;
			for (my $i = 0; $i < scalar @splitter; $i++){
				my $bin = ($i/2)+1;
				my $obs = $splitter[$i];
				$boot_obs_hash{$bin} = $splitter[$i];
				#print " i $i bin $bin value  $splitter[$i]\n";
				$hist_compl_obs[$bin] += $splitter[$i];
				
				$i++;
				my $exp = $splitter[$i];
				$boot_exp_hash{$bin} = $splitter[$i];
				$hist_compl_exp[$bin] += $splitter[$i];
				$cbo_minus_cbe_10bin_array[int($bin/10)] += $obs-$exp;
				#push @{$array_obs_minus_exp[$bin]}, $obs-$exp;
			}
			
			for (my $bin10 = 0; $bin10 < scalar @cbo_minus_cbe_10bin_array; $bin10++){
				push @{$array_cbo_minus_cbe[$bin10]}, $cbo_minus_cbe_10bin_array[$bin10];
			}
			## end of copypaste
			
			
			## 5.02 redundant
			#for (my $i = 0; $i < scalar @splitter; $i++){
			#	my $bin = ($i/2)+1;
			#	$boot_obs_hash{$bin} = $splitter[$i];
			#	$i++;
			#	$boot_exp_hash{$bin} = $splitter[$i];
			#}
			
			my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash);
			my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash);
			my $boot_obs_mean = hist_mean_for_hash(\%boot_obs_hash);
			my $boot_exp_mean = hist_mean_for_hash(\%boot_exp_hash);			
			$complement_boot_medians[$iteration][0] = $boot_obs_median;
			$complement_boot_medians[$iteration][1] = $boot_exp_median;
			$complement_boot_means[$iteration][0] = $boot_obs_mean;
			$complement_boot_means[$iteration][1] = $boot_exp_mean;
			print FILE "\n boot obs median: $boot_obs_median boot exp median $boot_exp_median  boot obs mean: $boot_obs_mean boot exp mean $boot_exp_mean\n";
			
			$iteration++;
		}
		close CSVFILE;
		
		
		
		my @array_diffdiff;
		if (scalar @array_gbo_minus_gbe  != scalar @array_cbo_minus_cbe){
			print "bintest failed! number of bin in group array is ".scalar @array_gbo_minus_gbe.", in complement array  ".scalar @array_cbo_minus_cbe."\n";
			print "Don't worry though, their equality is not required\n";
		}
		my $maxbin  = max(scalar @array_gbo_minus_gbe, scalar @array_cbo_minus_cbe);
		for (my $j = 0; $j < $maxbin; $j++){
			if (scalar @{$array_gbo_minus_gbe[$j]}  != scalar @{$array_cbo_minus_cbe[$j]}){
				print "itertest failed! number of iterations in group array for $j bin is ".scalar @{$array_gbo_minus_gbe[$j]}.", in complement array  ".scalar @{$array_cbo_minus_cbe[$j]}."\n";
				print "Now you can start worrying\n";
			}
			foreach my $iteration(@{$array_gbo_minus_gbe[$j]}){
				push @{$array_diffdiff[$j]}, $array_gbo_minus_gbe[$j][$iteration] - $array_cbo_minus_cbe[$j][$iteration];
			}
		}
		
		
		print FILE "bin mean_group_boot_obs mean_group_boot_exp ";
		print FILE "mean_compl_boot_obs mean_compl_boot_exp ";
		print FILE " diff_gbo-gbe-cbo+cbe_percentile_5 diff_gbo-gbe-cbo+cbe_percentile_95";
		print FILE " diff_gbo-gbe_percentile_5 diff_gbo-gbe_percentile_95";
		print FILE " diff_cbo-cbe_percentile_5 diff_cbo-cbe_percentile_95\n";
		my $maxbin  = max(scalar @hist_obs, scalar @hist_compl_obs);
		$maxbin = max ($maxbin, scalar @hist_exp);
		$maxbin = max ($maxbin, scalar @hist_compl_exp);
		
		for (my $j = 0; $j < $maxbin; $j++){ #foreach bin
			my $mean_group_obs = $hist_obs[$j]/$iteration;
			my $mean_group_exp = $hist_exp[$j]/$iteration;
			my $mean_compl_obs = $hist_compl_obs[$j]/$iteration;
			my $mean_compl_exp = $hist_compl_exp[$j]/$iteration;
			my $stat_gbo_minus_gbe = Statistics::Descriptive::Full->new();
			$stat_gbo_minus_gbe->add_data(\@{$array_gbo_minus_gbe[$j]});
			my $stat_cbo_minus_cbe = Statistics::Descriptive::Full->new();
			$stat_cbo_minus_cbe->add_data(\@{$array_cbo_minus_cbe[$j]});	
			my $stat_diffdiff = Statistics::Descriptive::Full->new();
			$stat_diffdiff->add_data(\@{$array_diffdiff[$j]});			
			print FILE "$j $mean_group_obs $mean_group_exp ";
			print FILE "$mean_compl_obs $mean_compl_exp ";
			print FILE $stat_diffdiff->percentile(5)." ".$stat_diffdiff->percentile(95);
			print FILE $stat_gbo_minus_gbe->percentile(5)." ".$stat_gbo_minus_gbe->percentile(95);
			print FILE $stat_cbo_minus_cbe->percentile(5)." ".$stat_cbo_minus_cbe->percentile(95)."\n";
		}
		
		
		
		my $pval_env_enrichment;
		my $pval_epi_enrichment;
		my $pval_env_depletion;
		my $pval_epi_depletion;		
		my $pval_epi;
		my $pval_env;
		my $pval_env_enrichment_for_mean;
		my $pval_epi_enrichment_for_mean;
		my $pval_env_depletion_for_mean;
		my $pval_epi_depletion_for_mean;		
		my $pval_epi_for_mean;
		my $pval_env_for_mean;
		
		for (my $i = 0; $i < scalar @group_boot_medians; $i++){
		## medians
			if (($group_boot_medians[$i][0] - $group_boot_medians[$i][1]) - ($complement_boot_medians[$i][0] - $complement_boot_medians[$i][1])
				>= ($obs_median - $exp_median) - ($complement_obs_median - $complement_exp_median)){
				$pval_env_enrichment += 1;
			}
			if (($group_boot_medians[$i][0] - $group_boot_medians[$i][1]) - ($complement_boot_medians[$i][0] - $complement_boot_medians[$i][1])
				<= ($obs_median - $exp_median) - ($complement_obs_median - $complement_exp_median)){
				$pval_env_depletion += 1;
			}
			if (-($group_boot_medians[$i][0] - $group_boot_medians[$i][1]) + ($complement_boot_medians[$i][0] - $complement_boot_medians[$i][1])
				>= -($obs_median - $exp_median) + ($complement_obs_median - $complement_exp_median)){
				$pval_epi_enrichment += 1;
			}
			if (-($group_boot_medians[$i][0] - $group_boot_medians[$i][1]) + ($complement_boot_medians[$i][0] - $complement_boot_medians[$i][1])
				<= -($obs_median - $exp_median) + ($complement_obs_median - $complement_exp_median)){
				$pval_epi_depletion += 1;
			}
			
			if ($group_boot_medians[$i][0] - $group_boot_medians[$i][1] >= $obs_median - $exp_median){
				$pval_env += 1;
			}
			if ($group_boot_medians[$i][0] - $group_boot_medians[$i][1] <= $obs_median - $exp_median){
				$pval_epi += 1;
			}
			
			## means
			
			if (($group_boot_means[$i][0] - $group_boot_means[$i][1]) - ($complement_boot_means[$i][0] - $complement_boot_means[$i][1])
				>= ($obs_mean - $exp_mean) - ($complement_obs_mean - $complement_exp_mean)){
				$pval_env_enrichment_for_mean += 1;
			}
			if (($group_boot_means[$i][0] - $group_boot_means[$i][1]) - ($complement_boot_means[$i][0] - $complement_boot_means[$i][1])
				<= ($obs_mean - $exp_mean) - ($complement_obs_mean - $complement_exp_mean)){
				$pval_env_depletion_for_mean += 1;
			}
			if (-($group_boot_means[$i][0] - $group_boot_means[$i][1]) + ($complement_boot_means[$i][0] - $complement_boot_means[$i][1])
				>= -($obs_mean - $exp_mean) + ($complement_obs_mean - $complement_exp_mean)){
				$pval_epi_enrichment_for_mean += 1;
			}
			if (-($group_boot_means[$i][0] - $group_boot_means[$i][1]) + ($complement_boot_means[$i][0] - $complement_boot_means[$i][1])
				<= -($obs_mean - $exp_mean) + ($complement_obs_mean - $complement_exp_mean)){
				$pval_epi_depletion_for_mean += 1;
			}
			
			if ($group_boot_means[$i][0] - $group_boot_means[$i][1] >= $obs_mean - $exp_mean){
				$pval_env_for_mean += 1;
			}
			if ($group_boot_means[$i][0] - $group_boot_means[$i][1] <= $obs_mean - $exp_mean){
				$pval_epi_for_mean += 1;
			}
		}
		
		print FILE "- pvalue_epistasis_enrichment pvalue_environment_enrichment pvalue_epistasis pvalue_environment\n";
		print FILE "median_stat ".($pval_epi_enrichment/$iteration)." ".($pval_env_enrichment/$iteration)." ".($pval_epi/$iteration)." ".($pval_env/$iteration)."\n";
		print FILE "mean_stat ".($pval_epi_enrichment_for_mean/$iteration)." ".($pval_env_enrichment_for_mean/$iteration)." ".($pval_epi_for_mean/$iteration)." ".($pval_env_for_mean/$iteration)."\n";

		close FILE;	
	
	
	}
	
}
close COUNTER;
}



sub no_check{
	return 1;
}


# prints protein_for_LRT files
sub print_data_for_LRT {
	$self = shift;
	my $dir = File::Spec->catdir(getcwd(),"likelihood", $self->{static_state});
	my $filename = File::Spec->catfile($dir, ($self->{static_protein})."_for_LRT.csv");
	open FILE, ">$filename" or die "Cannot create $filename";
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my @group = (1..565);
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($root);
	$self->visitor_coat ($root, \@array,\&lrt_visitor,\&no_check,\@args,0);
	print FILE "site,ancestor_node,t_branch_start,t_branch_end,event_indicator\n";
	foreach my $ind (@group){
		foreach my $ancnode(@{$self->{static_nodes_with_sub}{$ind}}){
			if(ref($ancnode) eq "REF"){
				$ancnode = ${$ancnode};
			}
			my $ancnodename = $ancnode->get_name();
			foreach my $node ( keys %{$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}}){
				print FILE $ind.",".$ancnodename.",".$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[1].",".
				$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[2].",";
				my $event = 0;
				if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]) {$event = 1};
				print FILE "$event\n";
			}
		}
	}	
	close FILE;
}




## Since 5.02
sub depth_groups_entrenchment_optimized_selection_alldepths {
	my $self = shift;
	my $step = $_[0];
	my $ancestral_nodes = $_[2];
	my $overwrite = $_[3];
	
	my $filename = $self->{static_output_base}.$self->{static_protein}."_for_enrichment";
	open FILE, ">>$filename";
	my $root = $self->{static_tree}-> get_root;
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
	my @args = ( $step, $root);
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0,$overwrite);

	foreach my $ind (@group){
	#print "here my ind $ind\n";
		foreach my $node(@{$self->{static_nodes_with_sub}{$ind}}){
			
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
				foreach my $bin (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
					$totmut += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$totlen += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
				}	
			#	print " totmut $totmut, totlen $totlen\n";
				next;
			}
		#	print "$site_node is still here\n";
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
	#print "maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
			if ($self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > 50){
				
				my %subtract_hash;
				
				
				if ($self->{static_subtract_tallest}){
					my $tallest_tip = ${$self->{static_hash_of_nodes}{$self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node; 
					my $path_length = $tallest_tip->get_branch_length; 
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last; # even if ancestors (from get_ancestors) contain the node itself, it won't make any difference
						}
						my $depth = $self->{static_distance_hash}{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				
				foreach my $bin (sort {$a <=> $b} keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
					#print "bin $bin adding ".$static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[0]."\n";
					$total_muts += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0){
					#print "total muts $total_muts \n";
					print FILE "site $ind node ".$node->get_name()." maxdepth ".$self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}."\n";
					my $check_local_lengths_sum;
					my $check_total_obs;
					my $check_total_exp;
					foreach my $bin (sort {$a <=> $b} (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}})){
						
							my $local_length = $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$check_local_lengths_sum += $local_length;
							$hist{$site_node}{$bin}[0] += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
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



# 21.12 Differs from depth_groups_entrenchment_optimized_selector_alldepths in that it keeps both site and node in obs_hash
# (yes, that's one line) 

sub depth_groups_entrenchment_optimized_selector_alldepths_2 {
	my $self = shift;
	my $step = shift;
	my $restriction = shift;
	unless (defined $restriction) {$restriction = 0};

	my $root = $self ->{static_tree}-> get_root;
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
	my @args = ( $step, $root);
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
		foreach my $node(@{$self ->{static_nodes_with_sub}{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $site_node = $ind."_".$node->get_name();
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				
				my %subtract_hash;
				
				if ($self ->{static_subtract_tallest}){
					#print "just checking ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}."\n";
					my $tallest_tip = ${$self ->{static_hash_of_nodes}{$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = $self ->{static_distance_hash}{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				foreach my $bin (keys %{$self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($self ->{static_subtract_tallest} && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0 && $total_muts > 0){ # 21.10 added total_muts > 0
				#print "total muts $total_muts \n";
				print "site $ind node ".$node->get_name()." maxdepth ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}."\n";
					foreach my $bin (sort {$a <=> $b} (keys %{$self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}})){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						 if ($total_length > 0){ 
							my $local_length = $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$site_node}{$bin}[0] += $self ->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
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

#26.02 almost depth_groups_entrenchment_optimized_selector_alldepths_2, but prints only total nodecounts for each site_node
sub nodeselector {
	my $self = shift;
	my $step = $_[0];
	my $restriction = $_[1];

	my $root = $self->{static_tree}-> get_root;
	my @array;
	my %hist;
	
	my @group;
	if ($_[3]){
		@group = @{$_[3]};
	}
	else {
		@group = (1..565);
	}
	 my $name = $_[4];

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root);
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	my %sitecounts;
	my %mutcounts;
	foreach my $ind (@group){
		foreach my $node(@{$self->{static_nodes_with_sub}{$ind}}){
			if(ref($node) eq "REF"){
				$node = ${$node};
			}
			my $site_node = $ind."_".$node->get_name();
			my $total_muts;
			my $total_length;
			
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			if ($self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				
				my %subtract_hash;
				
				if ($self->{static_subtract_tallest}){
					#print "just checking ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth_node"}."\n";
					my $tallest_tip = ${$self ->{static_hash_of_nodes}{$self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth_node"}}};	
					my @path_to_tallest_tip = @{$tallest_tip->get_ancestors()};	
					my $index_of_ancestor_node;
					my $path_length = $tallest_tip->get_branch_length; # todo: check if ancestors contain the node itself
				
					for(my $n = 0; $n < scalar @path_to_tallest_tip; $n++){		
						if ($path_to_tallest_tip[$n]->get_name eq $node->get_name){
							$index_of_ancestor_node = $n;
							last;
						}
						my $depth = $self->{static_distance_hash}{$node->get_name()}{$path_to_tallest_tip[$n]->get_name()};
						$subtract_hash{bin($depth,$step)} += $path_to_tallest_tip[$n]->get_branch_length;
						$path_length += $path_to_tallest_tip[$n]->get_branch_length;
					}
				}
				
				foreach my $bin (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}}){
					$total_muts += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];
					$total_length += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
					if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
						$total_length -= $subtract_hash{$bin};
					}
				}	
				#print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$static_subtree_info{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0 && $total_muts > 0){ # 21.10 added total_muts > 0
				#print "total muts $total_muts \n";
				my $md = $self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"};
				if ($md > 50){
					$sitecounts{50}{$ind} = 1;
				}
				if ($md > 100){
					$sitecounts{100}{$ind} = 1;
				}
				if ($md > 150){
					$sitecounts{150}{$ind} = 1;
				}				
				print "site $ind node ".$node->get_name()." maxdepth ".$self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"};
				my %totalcounts; # 26.02 counting nodes in analysis	
					foreach my $bin (sort {$a <=> $b} (keys %{$self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}})){
						#if ($total_length > 0 && $static_subtree_info{$node->get_name()}{$ind}{"hash"}{$bin}[1] > 0){ #there are some internal nodes with 0-length terminal daughter branches
						 if ($total_length > 0){ 
							my $local_length = $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[1];
							if ($self->{static_subtract_tallest} && $subtract_hash{$bin}){
								$local_length -= $subtract_hash{$bin};
							}
							$hist{$site_node}{$bin}[0] += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0]; #observed
							$hist{$site_node}{$bin}[1] += $total_muts*$local_length/$total_length; #expected
							$totalcounts{$site_node} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							
							if ($md > 50){
								$mutcounts{50} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							}	
							if ($md > 100){
								$mutcounts{100} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							}
							if ($md > 150){
								$mutcounts{150} += $self->{static_subtree_info}{$node->get_name()}{$ind}{"hash"}{$bin}[0];	# 26.02 counting nodes in analysis
							}							
						}
						if (!$hist{$site_node}{$bin}[0]){
							$hist{$site_node}{$bin}[0] += 0;
						}
						if (!$hist{$site_node}{$bin}[1]){
							$hist{$site_node}{$bin}[1] += 0;
						}

				}
					 print " ".$totalcounts{$site_node}."\n";
				}
			}
			
		}
	}	
	
		my $sites50 = scalar keys %{$sitecounts{50}};
		my $sites100 = scalar keys %{$sitecounts{100}};
		my $sites150 = scalar keys %{$sitecounts{150}};
		print $name." ".$sites50." ".$mutcounts{50}." ".$sites100." ".$mutcounts{100}." ".$sites150." ".$mutcounts{150}." "."\n";
	
	return %hist;	
}

sub my_median{
my @values = @{$_[0]};	
my $median;
my $mid = int @values/2;
my @sorted_values = sort @values;
if (@values % 2) {
    $median = $sorted_values[ $mid ];
} else {
    $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
} 
return $median;
}




#takes an array of probabilities for 0,1,2...
sub hist_median{
	my @hist = @{$_[0]};
	my $summ = sum (@hist);
	my $head = 0;
	my $interval = 0;
	my $median = 0;
	
	while ($head < $summ/2){
		$head += $hist[$interval];
		$median = $interval;
		$interval++;
	}
	
	if ($head == $summ/2){
		$median += 0.5;
	}
#print_hist(\@hist);
	return $median;
}

sub hist_mean {
	my @hist = @{$_[0]};
	my $summ = sum (@hist);
	my $integer;
	for(my $i = 0; $i <scalar @hist; $i++){
		$integer += $i*$hist[$i];
	}
	if ($summ > 0){
		return $integer/$summ;
	}
	else {
		return "NaN";
	}
}

#takes a hash of probabilities for 0,1,2...
sub hist_median_for_hash{
	my @hist = hist_to_array($_[0]);
	return hist_median(\@hist);
}

#takes a hash of probabilities for 0,1,2...
sub hist_mean_for_hash {
	my @hist = hist_to_array($_[0]);
	return hist_mean(\@hist);
}

# takes hash of probabilities for 0,1,2... and returns an array of ordered values
sub hist_to_array {
	my %prehist =  %{$_[0]};
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i};
		}
		else {
			push @hist, 0;
		}
	}
	return @hist;
}

sub hist_median_for_hash_arr{
	my %prehist =  %{$_[0]};
	my $number = $_[1];
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i}[$number];
		}
		else {
			push @hist, 0;
		}
	}
	return hist_median(\@hist);
}


sub hist_mean_for_hash_arr{
	my %prehist =  %{$_[0]};
	my $number = $_[1];
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	for (my $i = 1; $i <= $sorted_keys[-1]; $i++){
		if ($prehist{$i}){
			push @hist, $prehist{$i}[$number];
		}
		else {
			push @hist, 0;
		}
	}

	return hist_mean(\@hist);
}


## for hist->interval->site_index
sub hist_median_group {
	my @pre_hist = @{$_[0]};
	my @group = @{$_[1]};
	
	my @hist;
	foreach my $ind(@group){
		for (my $interval = 0; $interval < scalar @pre_hist; $interval++){
			$hist[$interval] += $pre_hist[$interval]->[$ind];
		}
	}
	
	return hist_median(\@hist);
}

sub print_hist {
	my @hist = @{$_[0]};

	my $counter = 0;
	print "\n";
	for (my $interval = 0; $interval < scalar @hist; $interval++){
			print $hist[$interval]."\t";
			$counter+=$hist[$interval];
		}
	print "\n";

}




sub median_difference{
	my @a1 = @{$_[0]};
	my @a2 = @{$_[1]};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a1);
	my $median1 = $stat->median();

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@a2);
	my $median2 = $stat->median();
	my $diff = $median1-$median2;
	print "median1 $median1 median2 $median2 diff $diff\n";
	return $diff;
}



1;
