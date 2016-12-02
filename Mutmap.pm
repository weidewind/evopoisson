#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

## 

package Mutmap;

use strict;
use Bio::Phylo::IO;


use Bio::Tools::CodonTable;
use TreeUtils::Phylo::FigTree;
use Bio::SeqIO;
use Data::Dumper;
use Bit::Vector;
use Try::Tiny;
use List::Util qw(sum min);
use Const::Fast;
use Switch;
use List::Util qw/shuffle/; 
use Statistics::Basic qw(:all);
use Statistics::TTest;
use Statistics::Descriptive;
use Storable qw(store retrieve lock_retrieve);
use IPC::System::Simple qw(capture);
use Class::Struct;
use IO::Handle;
use Cwd qw(abs_path cwd getcwd);
use File::Spec;
use Clone 'clone';
use Sub::Identify ':all';
use File::Path qw(make_path remove_tree);
use autodie;

#use DnaUtilities::observation_vector qw(make_observation_vector shuffle_obsv);
use observation_vector qw(make_observation_vector shuffle_obsv);
#use DnaUtilities::compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing);
use compare qw(nsyn_substitutions syn_substitutions nsyn_substitutions_codons is_neighbour_changing get_synmuts);
#use TreeUtils::Phylo::PhyloUtils qw(remove_zero_branches);
use PhyloUtils qw(remove_zero_branches);
use Groups;
use Memusage;
use Codeversion;

$| = 1;

	
	

	sub set_node_distance {
		my $self = shift;
		$self->{static_distance_hash}{$_[0]}->{$_[1]} = $_[2];	
	}
	
	sub set_alignment_length {
		my $self = shift;
		$self->{static_alignment_length} = $_[0]; 
	}
	
	sub has_node_distance {
		my $self = shift;
		if (!defined $self->{static_distance_hash}{$_[0]}->{$_[1]}){
			return 0;
		}
		else {
			return 1;
		}
	}
	
	sub get_node_distance {
		my $self = shift;
		return $self->{static_distance_hash}{$_[0]}->{$_[1]};
	}
	
	
	struct Mutation => {
		site_index => '$',
		node => '$',
	};
	
	
	sub maxpath_tag{
		my $subtract_maxpath = $_[0];
		my $tag;
		if (defined $subtract_maxpath){
			if ($subtract_maxpath eq "y" || $subtract_maxpath eq "yes" || $subtract_maxpath == 1 ){
				$tag = "maxpath_subtracted";
			}
			elsif ($subtract_maxpath eq "n" || $subtract_maxpath eq "no" || $subtract_maxpath == 0 ) {
				$tag = "maxpath_not_subtracted";
			}
			else {die "Invalid subtract_maxpath: $subtract_maxpath";}
		}
		else {$tag = '';}
		
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

	sub temp_tag {
			return "unreadable";
	}
	
	
	sub neighbour_tag {
		my $no_neighbour_changing = shift;
		my $tag;
		if (defined $no_neighbour_changing){
			if ($no_neighbour_changing eq "y" || $no_neighbour_changing eq "yes" || $no_neighbour_changing == 1 ){
				$tag = "no_neighbour_changing";
			}
			elsif ($no_neighbour_changing eq "n" || $no_neighbour_changing eq "no" || $no_neighbour_changing == 0 ) {
				$tag = "with_neighbour_changing";
			}
			else {die "Invalid no_neighbour_changing: $no_neighbour_changing";}
		}
		else {$tag = '';}
		
		return $tag;
	}
	
	sub leaves_tag {
		my $no_leaves = shift;
		my $tag;
		if (defined $no_leaves){
			if ($no_leaves eq "y" || $no_leaves eq "yes" || $no_leaves == 1 ){
				$tag = "no_leaves";
			}
			elsif ($no_leaves eq "n" || $no_leaves eq "no" || $no_leaves == 0 ) {
				$tag = "with_leaves";
			}
			else {die "Invalid no_leaves: $no_leaves";}
		}
		else {$tag = '';}
		
		return $tag;
	}

	sub printFooter {
		my $self = shift;
		my $outputStream = shift;
		print $outputStream "## protein ".$self->{static_protein}."\n";
		print $outputStream "## subtract_tallest ".$self->{static_subtract_tallest}."\n";
		print $outputStream "## state ".$self->{static_state}."\n";
		print $outputStream "## omit neighbour-changing mutations (only for 'reversals', ancestor n-ch muts are not skipped. Only valid for syn state)? 1 if true ".$self->{static_no_neighbour_changing}."\n";
		print $outputStream "## omit mutations on terminal branches? 1 if true ".$self->{static_no_leaves}."\n";
		print $outputStream "## output_base ".$self->{static_output_base}."\n";
		if ($self->{realdata}){
			print $outputStream "## realdata restriction ".get_realdata_restriction($self->{realdata})."\n";
		}
		print $outputStream "## code version hash ".Codeversion->get_version()."\n";
	}
	
	sub pathFinder {
		my $args = shift;	
		my $output_base = File::Spec->catdir(getcwd(), "output", $args->{bigdatatag}, $args->{bigtag}, state_tag($args->{state}), maxpath_tag($args->{subtract_tallest}), neighbour_tag($args->{no_neighbour_changing}), leaves_tag($args->{no_leaves})); 
		return $output_base;

	}
	
	sub createCodeversionFile {
		my $self = shift;
		my $script_name = shift;
		my $version_file = File::Spec->catfile($self->{static_output_base}, $script_name."_codeversion");
		open FILE, ">$version_file" or die "Cannot open $version_file: $!\n";
		print FILE (Codeversion::get_version());
		close FILE;
	}
	
	sub realdata_exists {
		my $args = shift;
		my $output_base = pathFinder($args);
		my $realdatapath = File::Spec->catfile($output_base, $args->{protein}."_".state_tag($args->{state})."_realdata");
		print "Checking if $realdatapath exists..";
		if (-f $realdatapath) {print "yes!\n";}
		else {print "no.\n";}
		return (-f $realdatapath);
	}
	
	# not to be confused with get_realdata_restriction, which needs realdata as an argument
	sub check_realdata_restriction{
		my $args = shift;
		my $output_base = pathFinder($args);
		my $realdatapath = File::Spec->catfile($output_base, $args->{protein}."_".state_tag($args->{state})."_realdata");
		my $realdata = lock_retrieve ($realdatapath) or die "Cannot retrieve ".$realdatapath;
		return get_realdata_restriction($realdata);
	}
	
	sub dataFinder {
		my $args = shift;	
		my $input_base = File::Spec->catdir(getcwd(), "data", $args->{bigdatatag});
		return $input_base;
	}

	sub set_tag {
		my $self = shift;
		my $tag = shift;
		my $output_subfolder = File::Spec->catdir($self->{static_output_base}, $tag);
		$self->{static_output_subfolder} = $output_subfolder;
		make_path($output_subfolder);
		make_path(File::Spec->catdir($output_subfolder, temp_tag()));
	}
	
	sub new {
		my ($class, $args) = @_;	
		#my $output_base = File::Spec->catdir(getcwd(), "output", $args->{bigdatatag}, $args->{bigtag}, state_tag($args->{state}), maxpath_tag($args->{subtract_tallest})); 
		my $output_base = pathFinder ($args);
		my $input_base = dataFinder ($args);
		my $treefile = File::Spec->catfile($input_base, $args->{protein}.".l.r.newick");
		my $static_tree = parse_tree($treefile)  or die "No tree at $treefile";
		my $self;
		make_path($output_base);
		make_path(File::Spec->catdir($output_base, temp_tag()));
		
		if ($args->{fromfile}){
			my $realdatapath;
			my $realdataname = $args->{protein}."_".state_tag($args->{state})."_realdata";
			if ($args->{fake}){
				 $realdatapath = File::Spec->catfile($output_base, $args->{tag}, $realdataname);
			}
			else { $realdatapath = File::Spec->catfile($output_base, $realdataname); }
			print "Creating mutmap from realdata $realdatapath\n";
			my $realdata = lock_retrieve ($realdatapath) or die "Cannot retrieve ".$realdatapath;
			if ($args->{fake}){
				unlink $realdatapath or warn "Could not unlink $realdatapath: $!";;
			}
			
			#foreach my $ind(1..300){
			#if ($realdata->{static_nodes_with_sub}{$ind}){
			#my $debugnum = scalar @{$realdata->{static_nodes_with_sub}{$ind}};
			#print "Early news from nnew: numnodes for $ind is $debugnum\n";
			#}
			#if ($realdata->{static_nodes_with_sub}{$ind} && scalar @{$realdata->{static_nodes_with_sub}{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
			#foreach my $node(@{$realdata->{static_nodes_with_sub}{$ind}}){
		#		print "Early News from nnew: nnode_name ".$$node->get_name()."\n";
		#	}
		#	}
		#	}
			
			
			$self = { 
				static_output_base => $output_base,
				static_input_base => $input_base,
				static_protein => $args->{protein},
				static_subtract_tallest => $args->{subtract_tallest},
				static_tree => $static_tree,
				static_treefile => $treefile,
				static_state => $args->{state},
				static_no_neighbour_changing =>$realdata->{no_neighbour_changing}, 
				static_no_leaves =>$realdata->{no_leaves},
				static_alignment_length => $realdata->{alignment_length}, 
				static_hash_of_nodes => $realdata->{hash_of_nodes}, 
				static_distance_hash => $realdata->{distance_hash},
				static_subs_on_node => $realdata->{static_subs_on_node}, # we never use these two when we produce new mutmappers from file (they are taken from observaton_vectors)
				obs_vectors => $realdata->{obs_vectors}, #added on 17.11.2016
				static_nodes_with_sub => $realdata->{static_nodes_with_sub}, #
				static_background_subs_on_node => $realdata->{bkg_subs_on_node},
				static_background_nodes_with_sub => $realdata->{bkg_nodes_with_sub},
				realdata => $realdata,
			};
			
			
		#	foreach my $ind(1..300){
		#	if ($self->{static_nodes_with_sub}{$ind}){
		#	my $debugnum = scalar @{$self->{static_nodes_with_sub}{$ind}};
		#	print "news from nnew: numnodes for $ind is $debugnum\n";
		#	}
		#	if ($self->{static_nodes_with_sub}{$ind} && scalar @{$self->{static_nodes_with_sub}{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
		##	foreach my $node(@{$self->{static_nodes_with_sub}{$ind}}){
		#		print "News from nnew: nnode_name ".$$node->get_name()."\n";
		#	}
		#	}
		#	}
			
		}
		else {
			my @arr = parse_fasta(File::Spec->catfile($input_base, $args->{protein}.".all.fa"));
			my %fasta = %{$arr[0]};
			my $alignment_length = $arr[1];
			my $static_protein  = $args->{protein};
			my %static_fasta = %fasta;
			my %static_hash_of_nodes;	
			my @nodes = $static_tree -> get_nodes;
			
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
				die "only syn or nsyn can be used as the second argument; unknown ".$args->{state}." was used instead";
			}

			$self = {
				static_output_base => $output_base,
				static_input_base => $input_base,
				static_protein => $args->{protein},
				static_alignment_length => $alignment_length, 
				static_subtract_tallest => $args->{subtract_tallest},
				static_no_neighbour_changing =>  $args->{no_neighbour_changing},
				static_no_leaves =>$args->{no_leaves},
				static_tree => $static_tree,
				static_treefile => $treefile,
				static_fasta => { %static_fasta },
				static_state  => $args->{state},
				static_hash_of_nodes => { %static_hash_of_nodes },
				static_subs_on_node => $mutmaps[0],
				static_nodes_with_sub => $mutmaps[1],
				static_background_subs_on_node => $bkg_mutmaps[0],
				static_background_nodes_with_sub => $bkg_mutmaps[1],
			};
			foreach my $node(@nodes){
				#if ($node->is_root()) {next;}
				my $name = $node ->get_name();
				$self ->{static_hash_of_nodes}{$name} = \$node;
			}
		}	
		
		bless $self, $class;
		$self->set_tag($args->{tag});
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
		my %nsyn = compare::nsyn_substitutions_codons($nodeseqs{$node->get_ancestors()->[0]->get_name()},
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
		my %nsyn = compare::syn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
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

sub mylength {
	my $self = shift;	
	my $length;
	#if ($self->{static_state} eq "nsyn"){ # 12.10 syn positions also correspond to codons
		$length = ($self->{static_alignment_length})/3;
	#	print "debugging mylength is $length\n";
	#}
	#elsif ($self->{static_state} eq "syn") {
	#	$length = $self->{static_alignment_length};
	#	print "debugging mylength is $length\n";
	#}
	#print "debugging returning mylength $length\n";
	return $length;
}
# sets static_sorted_nodnames and static_sorted_sites, retruns incidence_hash
sub incidence_matrix {
	my $self = shift;	
	my %matrix;
	my $length = $self->mylength();
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
		if ($self->{static_nodes_with_sub}{$ind}){
			my $debugnum = scalar @{$self->{static_nodes_with_sub}{$ind}};
			print "news from incidence: numnodes for $ind is $debugnum\n";
		}
		if ($self->{static_nodes_with_sub}{$ind} && scalar @{$self->{static_nodes_with_sub}{$ind}} > 0){# possibility to have at least 2 mutations after an ancestral one
			my %site_incidence = %empty_nodes_hash;
			push @sorted_sites, $ind;
			#print " added $ind to sorted sites\n";
			foreach my $node(@{$self->{static_nodes_with_sub}{$ind}}){
				print "News from incidence_matrix: nnode_name ".$$node->get_name()."\n";
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
	unless ($self->{static_sorted_nodnames} && $self->{static_sorted_sites}) {die "There is no static_sorted_nodnames (or static_sorted_sites) in this mutmap. Where did you take incidence_hash from?\n";}
	
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
    	my $self = shift;
    	my ( $node, $other_node ) = @_;
    	if  ($self->has_node_distance($node, $other_node)){
    		return $self->get_node_distance($node, $other_node);
    	}
    	else {
    		## calc_true instead of calc_my since 02 06 2015
    		my $dist = calc_true_patristic_distance($node, $other_node);
    		$self->set_node_distance($node, $other_node, $dist);
    		return $dist;
    	}
    	
    }


# prints tree with all mutations in the subtree of specified mutation (site, node). 
# If there is no such mutation, warns and proceeds.

sub print_subtree_with_mutations {
	my $self = shift;
	my @muts = @{@_[0]};
	my $tag = $_[1];
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my $myCodonTable   = Bio::Tools::CodonTable->new();
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($root);
	$self->visitor_coat ($root, \@array,\&lrt_visitor,\&no_check,\@args,0);
	my @output_files;
	for (my $i = 0; $i < scalar @muts; $i++){
		my ($ind, $ancnodename) = split(/_/, $muts[$i]);
		#my $ind = $muts[$i];
		#$i++;
		#my $ancnodename = $muts[$i];
		if (!exists $self->{static_subtree_info}{$ancnodename}{$ind}){
				warn "there is no mutation at $ind , $ancnodename";
		}
		
		my %sites;
		my %color;
		my $sub = ${$self -> {static_subs_on_node}{$ancnodename}}{$ind};
		$sites{$ancnodename} = $ancnodename."_".$ind."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
		$color{$ancnodename} = "-16776961";
		foreach my $node ( keys %{$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}}){
			if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]){
				my $sub = ${$self -> {static_subs_on_node}{$node}}{$ind};
				print ($sub."\n");
				$sites{$node} = $node."_".$ind."_".$sub->{"Substitution::ancestral_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::ancestral_allele"}).")->".
							  $sub->{"Substitution::derived_allele"}.
							  "(".$myCodonTable->translate($sub->{"Substitution::derived_allele"}).")";
				$color{$node} = "-16776961";
			}
			else {
				$color{$node} = "-16776961";
			
				#print $file $ind.",".$ancnodename.",".$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[1].",".
				#$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[2].",";
				#my $event = 0;
				#if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]) {$event = 1};
				#print $file "$event\n";
			}
		}
		my $file = $self -> {static_protein}."_sites_".$ind."_".$ancnodename.".tre";
		my $dir = File::Spec -> catdir($self -> {static_output_base}, "trees", $tag);
		make_path($dir);
		my $filepath = File::Spec -> catfile($dir, $file);
		open TREE, ">$filepath";
		print TREE "#NEXUS\n\nbegin trees;\n";
		print TREE "\ttree $ind = [&R] ";
		my $tree_name=tree2str($self -> {static_tree},sites => \%sites, color=>\%color);
		print TREE $tree_name;
		print TREE "\nend;\n";
		my $figblock = File::Spec -> catfile(getcwd(), "figtree_block");
		open BLOCK, "<$figblock" or die "Cannot open figtree_block: $!\n";
		while (<BLOCK>){
			print TREE $_;
		}
		close BLOCK;
		close TREE;
		push @output_files, $filepath;
	}
	return @output_files;
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
	
	my $file = $self -> {static_protein}."_sites_".$site.".tre";
	my $dir = File::Spec -> catdir($self -> {static_output_base}, "trees");
	make_path($dir);
	my $filepath = File::Spec -> catfile($dir, $file);
	open TREE, ">$filepath";
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





sub myclone {
	my $self = shift;
	my $clone = {
			static_output_base => $self->{static_output_base},
			static_protein => $self->{static_protein},
			static_tree =>  $self->{static_tree},
			static_fasta => $self->{static_fasta},
			static_state  => $self->{static_state},
			static_alignment_length  => $self->{static_alignment_length},
			static_hash_of_nodes => $self->{static_hash_of_nodes},
			static_distance_hash => $self->{realdata}{"distance_hash"},
			static_background_subs_on_node => $self->{static_background_subs_on_node },
			static_background_nodes_with_sub => $self->{static_background_nodes_with_sub},
			obs_vectors => clone($self->{realdata}{"obs_vectors"})  #the only structure which can (and will) be changed
	};
	
	bless $clone, ref $self; #ref $self returns class of object $self
	return $clone;
}



# outputs hash of hashes used for construction of observed_vector
sub get_hashes {
	my $self = shift;
	my %res_hash;
	print ("what must be a hashref is a ".ref($self->{static_nodes_with_sub})."\n");
	foreach my $ind(keys $self->{static_nodes_with_sub}){
		my %x_hash;
		my %y_hash;
		foreach my $node($self->{static_tree}->get_nodes){
			$x_hash{$node->get_name()} = $node->get_branch_length;
			if ($self->{static_subs_on_node}{$node->get_name()}{$ind}){
				$y_hash{$node->get_name()} = 1;			
			}
			else {
				$y_hash{$node->get_name()} = 0;
			}
		}
		$res_hash{$ind}{"x"} = \%x_hash;
		$res_hash{$ind}{"y"} = \%y_hash;
	}
	return %res_hash;	
}



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
#my $counter = 0;	
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
#$counter++;			
			push (@{$nodes_with_sub{$ind}}, \${$self->{static_hash_of_nodes}{$nodname}});
		#print "TEST1 ".${$static_hash_of_nodes{$nodname}}->get_name()."\n"; # часть имен исчезла, а часть - осталась ОО
		#print "TEST2 ".$nodname."\n";
		}
		}
	}
#print "There are $counter mutations in this clone\n";	
	return (\%subs_on_node, \%nodes_with_sub);
}



sub shuffle_mutator {
	my $self = shift;
	my %obs_vectors = $self->get_observation_vectors(); # created only once, reused afterwards
	my %shuffled_obs_vectors = shuffle_observation_vectors(\%obs_vectors);
	$self->{obs_vectors} = \%shuffled_obs_vectors;
	my @mock_mutmaps = $self->read_observation_vectors(\%shuffled_obs_vectors); 
	$self->{static_subs_on_node} = $mock_mutmaps[0];
	$self->{static_nodes_with_sub} = $mock_mutmaps[1];
	return $self; 
}

sub FDR_all {
	my $self = shift;
	my $number_of_fakes = shift;
	my $mock_mutmap = $self->myclone();
	for (my $i = 1; $i <= $number_of_fakes; $i++){
		$mock_mutmap->shuffle_mutator();
	}
}




# 5.11 for entrenchment_bootstrap_full_selection_vector
sub iterations_gulp {
	my $self = shift;
	my $iterations = shift;
	my $tag = shift;
	my $verbose = shift;
	my $memusage = shift;
	my $restriction = shift;
		
	#my $ancestor_nodes = $_[4];
	#my $obs_vector = $_[5];
	#my $norm = $_[6];
	if ($verbose){print "Extracting realdata..\n";}	
	my $realdata = $self->{realdata};
	my $maxbin = $realdata->{"maxbin"};
	my $step = $realdata->{"step"}; #bin size
	unless (defined $step) {die "Oh no, bin size in realdata is not defined. Won't proceed with simulations.\n";}
	my $ancestor_nodes = $realdata->{"ancestor_nodes"};
	#my $obs_vectors = $realdata->{"obs_vectors"};
	my $outdir;
	if ($self->{static_output_subfolder}){
		$outdir = $self->{static_output_subfolder};
	}
	else {
		$outdir = $self->{static_output_base};
	}
	if ($verbose){print "Cloning mutmap..\n";}
	my $mock_mutmap = $self->myclone(); # 25.07 create a new object for shuffling
	my @simulated_hists;
	
	for (my $i = 1; $i <= $iterations; $i++){
		#if ($verbose){print "Creating clone..\n";}
		#my $mock_mutmap = $self->myclone(); # 30.11 test
		if ($verbose){print "Shuffling clone..\n";}
		$mock_mutmap->shuffle_mutator(); # this method shuffles observation vectors and sets new $static_nodes.. and static_subs..
		my %hash;
		# >new iteration string and all the corresponding data  are printed inside this sub:
		my %prehash = $mock_mutmap->depth_groups_entrenchment_optimized_selection_alldepths($step,$restriction,$ancestor_nodes, "overwrite", $tag, $verbose); #step (bin size), restriction, ancestor_nodes, should I overwrite static hash?

		foreach my $bin(1..$maxbin){
				foreach my $site_node(keys %prehash){
					$hash{$bin}[1] += $prehash{$site_node}{$bin}[1];
					$hash{$bin}[0] += $prehash{$site_node}{$bin}[0];				
			}
		}
		# 21.12 we do not need any norm here since we normalize values in concat_and_divide - separately for different maxdepths
		push @simulated_hists, \%hash;
		#%static_ring_hash = (); # must be cleaned in visitor_coat
		#%static_subtree_info = (); # must be cleaned in visitor_coat
	}
	if ($memusage){
		my $locker = Memusage->get_locker($self);
		$locker->print_memusage();
	}
	# store \@simulated_hists, File::Spec->catfile($self->{static_output_base}, $self->{static_protein}."_gulpselector_vector_alldepths_stored_".$tag); # was used because I was afraid of loosing a large amount of time because of some mistake
	# my $arref = retrieve(File::Spec->catfile($self->{static_output_base}, $self->{static_protein}."_gulpselector_vector_alldepths_stored_".$tag));
	my $arref = \@simulated_hists;
	my $csvfile = File::Spec->catfile($outdir, temp_tag(), $self->{static_protein}."_gulpselector_vector_alldepths_".$tag.".csv");
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



sub get_obshash {
	my $realdata = shift;
	my $restriction = shift;
	my $rr = get_realdata_restriction($realdata);
	print "realdata restriction is $rr and we need $restriction\n"; 
	unless(defined $rr && $rr <= $restriction ){
			die "realdata restriction is undefined or is greater than get_obshash restriction: ".$rr." > $restriction \n";
	}
	my $obs_hash = $realdata->{"obs_hash".$rr};
	return $obs_hash;
}

# not to be confused with check_realdata_restriction, which gets constructor arguments as an argument
sub get_realdata_restriction {
	my $realdata = shift;
	#my @obshash_restriction = map { /^obs_hash(.*)/ ? $1 : () } (keys %{$realdata});
	#return $obshash_restriction[0];
	return $realdata->{restriction};
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
		@group = (1..$self->mylength());
	}
	my %group_hash;
	foreach my $ind(@group){
		$group_hash{$ind} = 1;
	}
	
#	$real_data = lock_retrieve("/export/home/popova/workspace/perlCoevolution/TreeUtils/Phylo/MutMap/".$prot."_realdata") or die "Cannot retrieve real_data";
	my $realdata = $self->{realdata};
	my $obs_hash = get_obshash($realdata, $restriction);
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


# 13.09.2016 compute_norm for one site_node (for single site poisson analysis)
sub compute_norm_single_site {
	my $self = shift;
	my $site_node = shift;

	my $realdata = $self->{realdata};
	my $obs_hash = get_obshash($realdata, 1000); # 1000 is supposed to be bigger than any restriction in realdata, so this will just silently give you obshash from realdata
	my $subtree_info = $realdata->{"subtree_info"};
	
	my $norm;
	foreach my $bin(keys %{$obs_hash->{$site_node}}){
		$norm += $obs_hash->{$site_node}->{$bin}->[0];
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

		#%static_ring_hash = ();
		#%static_subtree_info = ();
		$i++;
	}
}



#5.11 for entrenchment_bootstrap_full_selection_vector
# analyze real data, prepare bkg_mutmap, ancestor_nodes and obs_vector
#28.12 hash is pruned: we do not keep mutation info if its maxdepth<50
#1.08 hash is not necessarily pruned - you can set restriction to 0 to get complete data
sub prepare_real_data {
	my $self = shift;
	my $args = shift;
	my $restriction = $args->{restriction};
	my $step = $args->{step};
	my $fake = $args->{fake};
	unless(defined $step) { $step = 0.5; }
	unless(defined $restriction) { $restriction = 50; }
	my $prot = $self->{static_protein};
	$self -> set_distance_matrix();
	my %matrix = $self->incidence_matrix(); 
	$self -> print_incidence_matrix(\%matrix);
	my $debugnum = scalar keys %{$self ->{static_nodes_with_sub}};
	print "Very early News from prepare: static_nodes_with_sub contains $debugnum keys\n";
	# used depth_groups_entrenchment_optimized_selector_alldepths but changed it for depth_groups_entrenchment_optimized_selector_alldepths_2, because the latter
	# keeps in obs_hash info about site index as well as about node name
	my %full_obs_hash = $self -> depth_groups_entrenchment_optimized_selector_alldepths_2($step, $restriction); # bin size
	my %ancestor_nodes;
	foreach my $ancnode(keys %full_obs_hash){
	my @splitter = split(/_/, $ancnode);
		$ancestor_nodes{$splitter[-1]} = 1;
		print "Ancestor: ".$splitter[-1]."\n";
	}	
	my $restricted_norm;
	my %restricted_obs_hash;
	my $maxbin = 0;
	my $debugnum = scalar keys %full_obs_hash;
	print "Early news from prepare: there are $debugnum keys in full_obs_hash\n";
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
	my $debugnum = scalar keys %restricted_obs_hash;
	print "News from prepare: there are $debugnum keys in restricted_obs_hash\n";
	my %realdata = (
		"norm".$restriction => $restricted_norm,
		step => $step,
		restriction => $restriction,
		maxbin => $maxbin,
		ancestor_nodes => \%ancestor_nodes,
		obs_vectors => \%obs_vectors,
		bkg_subs_on_node => $self -> {static_background_subs_on_node},
		bkg_nodes_with_sub => $self -> {static_background_nodes_with_sub},
		distance_hash => $self -> {static_distance_hash},
		hash_of_nodes => $self -> {static_hash_of_nodes},
		subtree_info => $self -> {static_subtree_info},
		alignment_length => $self -> {static_alignment_length},
		no_neighbour_changing => $self -> {static_no_neighbour_changing},
		no_leaves => $self -> {static_no_leaves},
		"obs_hash".$restriction => \%restricted_obs_hash,
	);
	
	## added at 17.11.2016 for fake mutmaps (mutmap, produced from realdata, now can be used for printing it (these hashes are necessary for depth_.._2))
		$realdata{static_subs_on_node} = $self -> {static_subs_on_node}; # if it used for fake, it will 
		$realdata{static_nodes_with_sub} = $self -> {static_nodes_with_sub};
	##
	
	
	my $realdatapath;
	if ($fake){
		$realdatapath = $self->{static_output_subfolder};
	} 
	else {
		$realdatapath = $self->{static_output_base};
	}
	$realdatapath = File::Spec->catfile($realdatapath, $prot."_".state_tag($self->{static_state})."_realdata");
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
	#my @obshash_restriction = map { /^obs_hash(.*)/ ? $1 : () } (keys $realdata);
	#unless(defined $obshash_restriction[0] && $obshash_restriction[0] <= $restriction ){
	#		die "realdata restriction is bigger than select_ancestor_nodes restriction: ".$obshash_restriction[0]." > $restriction \n";
	#}
	#my $obs_hash = $realdata->{"obs_hash".$obshash_restriction[0]};
	my $obs_hash = get_obshash($realdata, $restriction);
	my $subtree_info = $realdata->{"subtree_info"};
	my %group_nodes;
	my $debugnum = scalar keys %{$obs_hash};
	print "there are $debugnum keys in obshash\n";
	foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = split(/_/, $site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
				if ($maxdepth > $restriction && $group_hash{$site}){
					$group_nodes{$node_name} = 1;
					print "group_node ".$node_name."\n";
				}
		}
		my $count = scalar keys %group_nodes;
	print "Total $count\n";
	return %group_nodes;
}	
	
	

	

sub count_iterations {
	my $self = shift;
	my $prot = $self->{static_protein};
	my $dir = $self->{static_output_base};
	my $dirname = File::Spec->catdir($dir, $prot); 
	make_path ($dirname);
	opendir(DH, $dirname);
	my @files = readdir(DH);
	closedir(DH);
	unless (scalar @files > 0){
		return 0;
	}
	my $counter = 0;
	foreach my $gulp_filename(@files){
		next if (-d $gulp_filename);
		my $fullpath = File::Spec -> catfile($dirname, $gulp_filename);
		open GULP, "<$fullpath" or die "Cannot open $fullpath";
		while(<GULP>){
			if ($_ =~ /^>.*/){$counter++;}
		}
	}
	return $counter;
}

sub iterations_maxtag {
	my $self = shift;
	my $prot = $self->{static_protein};
	my $dir = $self->{static_output_base};
	my $dirname = File::Spec->catdir($dir, $prot); 
	make_path ($dirname);
	opendir(DH, $dirname);
	my @files = readdir(DH);
	closedir(DH);
	unless (scalar @files > 0){
		return 0;
	}
	my @tags = 0;
	foreach my $gulp_filename(@files){
		next if (-d $gulp_filename);
		if ($gulp_filename =~ /.*_([0-9]+)$/){
			push @tags, $1;
		}
		else {
			print "Strange tag in iterations file $gulp_filename ! It might be an error, and it might not.\nWon't change my behavior because of some stupid tag, just wanted you to be aware of it.\n";

		}
	}
	my $maxtag = List::Util::max(@tags);
	return $maxtag;
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
	my $subdir = $self->{static_output_subfolder};
	if (! defined $subdir){
		$subdir = $dir;
	}
	my $nodecount_file = File::Spec->catfile($subdir, $prot."_nodecount");
#	open NODECOUNT, ">$nodecount_file" or die "Cannot create $nodecount_file";
	
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
	make_path ($dirname);
	opendir(DH, $dirname);
	my @files = grep { /.*_[0-9]+/ }readdir(DH);
	unless (scalar @files > 0){die "No simulation files found in folder $dirname\n";}
	closedir(DH);
	
	my %filehandles;

	foreach my $md(@maxdepths){
		foreach my $group_number(0.. scalar @groups - 1){
			local *FILE;
			my $csvfile =  File::Spec->catfile($subdir, temp_tag(),$prot."_gulpselector_vector_".$md."_".$group_names[$group_number].".csv");
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
				#print "summtest ok\n";
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
						#			print NODECOUNT "maxdepth $md group ".$group_names[$group_number]." node $node_name\n";
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
			#		print "going to print something\n";
					foreach my $bin(1..$maxbin){
						print $filehandle $hash{$md}[$group_number]{$bin}[0].",".$hash{$md}[$group_number]{$bin}[1].",";
					}
					print $filehandle "\n";
				
				
				}
			}
			

			
			$iteration_number++;
			if ($iteration_number%50 == 0){
				print $iteration_number."\n";
			}
		}
		
		close GULP;
		
	}
	
#	close NODECOUNT;
	
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
	

#13.09.2016 prints one file for one ancestor site_node (instead of group of nodes) 	
	
sub concat_and_divide_simult_single_sites {
	my $self = shift;
	my $prot = $self->{static_protein};
	my @maxdepths = @{$_[0]};
	#my @groups = @{$_[1]};
	#my @group_names = @{$_[2]};
	my $subtract_maxpath = $self->{static_subtract_tallest};
	my $dir = $self->{static_output_base};
	my $subdir = $self->{static_output_subfolder};
	if (! defined $subdir){
		$subdir = $dir;
	}
	my $realdata = $self->{realdata};
	my $obs_hash = get_obshash($realdata, List::Util::min(@maxdepths));
	my $subtree_info = $realdata->{"subtree_info"};
	
	my %hash;
	my $iteration_number = 1;
	my $maxbin = $realdata->{"maxbin"};
	
	my %norms; #previous: $norms{$md}[$group_number] = norm
	foreach my $md(@maxdepths){
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = split(/_/, $site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $md){
				$norms{$md}{$site_node} = $self->compute_norm_single_site($site_node);
			}
		}	
	}
	
	
	
	my $dirname = File::Spec->catdir($dir, $prot); 
	make_path ($dirname);
	opendir(DH, $dirname);
	my @files = grep { /.*_[0-9]+$/ }readdir(DH); 
	unless (scalar @files > 0){die "No simulation files found in folder $dirname\n";}
	closedir(DH);
	
	my %filehandles;

	foreach my $md(@maxdepths){
		foreach my $site_node(keys %{$obs_hash}){
			local *FILE;
			my $csvfile =  File::Spec->catfile($subdir, temp_tag(),$prot."_gulpselector_vector_".$md."_".$site_node.".csv");
			open FILE, ">$csvfile" or die "Cannot create $csvfile";
			FILE->autoflush(1);
			$filehandles{$md}{$site_node} = *FILE;
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
			my @str_array;
			my $str = <GULP>;
			
			my %sums;
			my %hash;
			
			my $test_obs_summ;
			my $test_exp_summ;
			
			while ($str =~ /^[^>]/){ 

				if ($str =~ /^site/){
					if ($test_obs_summ - $test_exp_summ < 0.0001 && -$test_obs_summ + $test_exp_summ < 0.0001){
					#	print "summtest ok\n";
					}
					else {
						print "summtest failed! $site $node_name obssum $test_obs_summ, expsum $test_exp_summ\n";
					}
					$test_obs_summ = 0;
					$test_exp_summ = 0;	
				
					my @str_array = split(/\s+/, $str);
					$site = $str_array[1]; # careful (never used afterwards, only for debugging)
					$node_name = $str_array[3];
					$max_depth = $str_array[5];
					
					
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
						foreach my $site_node(keys %{$obs_hash}){
							my ($mysite, $mynode) = split(/_/, $site_node);
							if ($mynode eq $node_name){ # 28.09.2016 
								if ($norms{$md}{$site_node}){ 
								#print "group number $group_number md $md node name $node_name\n";
									$sums{$md}{$site_node} += $str_array[1];
									$hash{$md}{$site_node}{$str_array[0]}[1] += $str_array[2];
									$hash{$md}{$site_node}{$str_array[0]}[0] += $str_array[1];
									#print $hash{$md}{$site_node}{$str_array[0]}[0]." obs integral\n";
								}
							}
						}
					}
				}

				
				$str = <GULP>;
				
			}
			## parsed all information about this iteration 

			# maxbins are different for every iteration. Find maximum and use it.

			$maxbin = max($maxbin, $str_array[0]);
#print "maxbin $maxbin, iteration number $iteration_number\n";	
#print "sum50 $sum50 sum100 $sum100 sum150 $sum150 norm 50 $norm50 norm 100 $norm100 norm 150 $norm150\n";		
			
			foreach my $md(@maxdepths){ 
				foreach my $site_node(keys %{$obs_hash}){
					#print "maxdepth $md group number $group_number \n";
					if ($sums{$md}{$site_node} == 0){
						foreach my $bin(1..$maxbin){
							$hash{$md}{$site_node}{$bin}[0] = "NA";
							$hash{$md}{$site_node}{$bin}[1] = "NA";
						}
					}
					else {
						foreach my $bin(1..$maxbin){
							#print "in hash: ".$hash{$md}[$group_number]{$bin}[0]."\n";
							#print "norm ".$norms{$md}[$group_number]."\n";
							#print "sum ".$sums{$md}[$group_number]."\n";
							$hash{$md}{$site_node}{$bin}[0] = $hash{$md}{$site_node}{$bin}[0]*$norms{$md}{$site_node}/$sums{$md}{$site_node};
							$hash{$md}{$site_node}{$bin}[1] = $hash{$md}{$site_node}{$bin}[1]*$norms{$md}{$site_node}/$sums{$md}{$site_node};
						}
					}
					
					my $filehandle = $filehandles{$md}{$site_node};
					#print "going to print something\n";
					foreach my $bin(1..$maxbin){
						print $filehandle $hash{$md}{$site_node}{$bin}[0].",".$hash{$md}{$site_node}{$bin}[1].",";
					}
					print $filehandle "\n";
				
				
				}
			}
			
			
			$iteration_number++;
			if ($iteration_number%50 == 0){
				print "iteration number ".$iteration_number."\n";
			}
		}
		
		close GULP;
		
	}
	
	
	foreach my $md(@maxdepths){
		foreach my $site_node(keys %{$obs_hash}){
					my $filehandle = $filehandles{$md}{$site_node};
					close $filehandle;
		}
		
	}
	

}	
	

# counter from count_pvalues
sub group_counter {
	my $self = $_[0];
	my $prot = $self -> {static_protein};
	#my $prot = $_[0];
	my @restriction_levels = @{$_[1]};
	my @groups = @{$_[2]};
	my @group_names = @{$_[3]};
	my $dir = $self -> {static_output_base};
	
		my $countfile = File::Spec->catfile($dir, $prot."_count");
	open COUNTER, ">$countfile" or die "Cannot create $countfile";
	COUNTER->autoflush(1);
	
	my $realdata =  $self -> {realdata};
	#my $realdata = get_real_data();
	
	my $maxbin = $realdata->{"maxbin"}; 
	my $step =  $realdata->{"step"};
	unless (defined $step) {die "Oh no, realdata bin size is not defined. Won't proceed with pvalues\n";}
	#print "before cycle\n";
	
	
	my $obs_hash = get_obshash($realdata, List::Util::min(@restriction_levels)); # if min $restriction is less than restriction in realdata, it will die
	my $subtree_info = $realdata->{"subtree_info"};
	for my $restriction(@restriction_levels){
		print "level $restriction\n";
		
		# only for @all@
		my $group_number = scalar @groups - 1;
		print " groups ".scalar @groups - 1;
		my %group_hash;
		print " size ".scalar @{$groups[$group_number]}."\n";
		foreach my $site(@{$groups[$group_number]}){
			#print "test-1 $site\n";
			$group_hash{$site} = 1;
		}
		my %obs_hash_restricted;
		my $norm_restricted;
		## copypaste from prepare: create restricted hash	
		# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
		my $mutcounter;
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = split(/_/, $site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $restriction && $group_hash{$site}){ #15.02 moved here
			foreach my $bin(keys %{$obs_hash->{$site_node}}){
				$maxbin = max($bin, $maxbin);
					$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
				}
		
			}
		#	print "MAXBIN $maxbin\n";
		}
		my $count = scalar keys %obs_hash_restricted;
		print COUNTER "$restriction all group $count muts $norm_restricted\n";
		print COUNTER "$restriction all group $count muts $norm_restricted\n";
			for (my $group_number = 0; $group_number < scalar @groups - 1; $group_number++){ 
	
			## group
			my %group_hash;
			
			foreach my $site(@{$groups[$group_number]}){
				$group_hash{$site} = 1;
			}
			my %obs_hash_restricted;
			my $norm_restricted;
			## copyaste from prepare: create restricted hash	
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
			#	print "MAXBIN $maxbin\n";
			}
			
			## end f copypaste	
			my $count = scalar keys %obs_hash_restricted;
			print  COUNTER "$restriction ".$group_names[$group_number]." group $count muts $norm_restricted\n";
			print  COUNTER "$restriction ".$group_names[$group_number]." group $count muts $norm_restricted\n";
			print COUNTER   "$restriction ".$group_names[$group_number]."\n";
			foreach my $sn(keys %obs_hash_restricted){
				print  COUNTER $sn."\n";
			}
			$group_number++;
			}
		
	}
	close COUNTER;
}	
	
	
# 19.12 after concat_and_divide	
sub count_pvalues{	
	my $self = $_[0];
	my $prot = $self -> {static_protein};
	#my $prot = $_[0];
	my @restriction_levels = @{$_[1]};
	my @groups = @{$_[2]};
	my @group_names = @{$_[3]};
	my $fake = $_[4];
	my $dir = $self -> {static_output_base};
	my $outdir = $self -> {static_output_subfolder};
	#my $dir = $_[5];
	
	my $countfile = File::Spec->catfile($outdir, $prot."_count");
	if ($fake) {
		open COUNTER, ">>$countfile" or die "Cannot create $countfile";
	}
	else {
		open COUNTER, ">$countfile" or die "Cannot create $countfile";
	}
	COUNTER->autoflush(1);
	
	my $realdata =  $self -> {realdata};
	#my $realdata = get_real_data();
	
	my $maxbin = $realdata->{"maxbin"}; 
	my $step =  $realdata->{"step"};
	unless (defined $step) {die "Oh no, realdata bin size is not defined. Won't proceed with pvalues\n";}
	#print "before cycle\n";
	
	
	my $obs_hash = get_obshash($realdata, List::Util::min(@restriction_levels)); # if min $restriction is less than restriction in realdata, it will die
	#my $obs_hash = $realdata->{"obs_hash50"}; 
	my $subtree_info = $realdata->{"subtree_info"};
	for my $restriction(@restriction_levels){
		print "level $restriction\n";
		
		# only for @all@
		my $group_number = scalar @groups - 1;
		print " groups ".scalar @groups - 1;
		my %group_hash;
		print " size ".scalar @{$groups[$group_number]}."\n";
		foreach my $site(@{$groups[$group_number]}){
			#print "test-1 $site\n";
			$group_hash{$site} = 1;
		}
		my %obs_hash_restricted;
		my $norm_restricted;
		## copypaste from prepare: create restricted hash	
		# Beware! obs_hash50 from real_data really contains restricted data (maxdepth>50)
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = split(/_/, $site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
		#	if (compare::is_neighbour_changing($self->$static_subs_on_node{$node_name}{$site}, 1) == 1) # fisk
			if ($maxdepth > $restriction && $group_hash{$site}){ #15.02 moved here
			foreach my $bin(keys %{$obs_hash->{$site_node}}){
				$maxbin = max($bin, $maxbin);
				#if ($maxdepth > $restriction && $group_hash{$site}){ #15.02 moved from here
					$norm_restricted += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
				}
		
			}
		#	print "MAXBIN $maxbin\n";
		}
		my $count = scalar keys %obs_hash_restricted;
		print COUNTER "Total for $restriction all $count\n";
		
			
		## end of copypaste	
			
		my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]);
		my $outputfile;
		if ($fake){
			open $outputfile, ">>$file" or die "Cannot create $file";
		}
		else {
			open $outputfile, ">$file" or die "Cannot create $file";
		}
				
		my %histhash;
		foreach my $site_node(keys %obs_hash_restricted){
			foreach my $bin (keys %{$obs_hash_restricted{$site_node}}){
				$histhash{$bin}[0] += $obs_hash_restricted{$site_node}{$bin}[0];
				$histhash{$bin}[1] += $obs_hash_restricted{$site_node}{$bin}[1];
			}
		}
		print $outputfile "bin\tobs\texp\n";
		my @sorted_bins = sort { $a <=> $b } keys %histhash;
		foreach my $bin (@sorted_bins){
			print $outputfile $bin."\t".$histhash{$bin}[0]."\t".$histhash{$bin}[1]."\n";
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
		my $obs_median = hist_median_for_hash(\%flat_obs_hash, $step);
		my $exp_median = hist_median_for_hash(\%flat_exp_hash, $step);
		my $obs_mean = hist_mean_for_hash(\%flat_obs_hash, $step); # 18.03 - added the same statistics based on histogram mean (instead of median)
		my $exp_mean = hist_mean_for_hash(\%flat_exp_hash, $step);
		
		print $outputfile "\n observed median: $obs_median\n";
		print $outputfile "\n poisson expected median: $exp_median\n";
		print $outputfile "\n observed mean: $obs_mean\n";
		print $outputfile "\n poisson expected mean: $exp_mean\n";
		
		my $pval_epi;
		my $pval_env;
		my $pval_epi_for_mean;
		my $pval_env_for_mean;
	
	
		#print "going  read input file\n";	
		
		if ($obs_mean ne "NaN"){
		
		my $csvfile = File::Spec->catfile($outdir, temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number].".csv");
		open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
		my $iteration = 0;
		my @hist_obs;
		my @hist_exp;
		my @array_obs_minus_exp;
		while(<CSVFILE>){
			my %boot_obs_hash;
			my %boot_exp_hash;
			my @splitter = split(/,/, $_);
			if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
				next;
			}
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
				push @{$array_obs_minus_exp[$bin10]},  $diff_10bin_array[$bin10]; # no itnumber needed
			}
			
			
			
			my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash, $step);
			my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash,$step);
			my $boot_obs_mean = hist_mean_for_hash(\%boot_obs_hash, $step);
			my $boot_exp_mean = hist_mean_for_hash(\%boot_exp_hash, $step);
			print $outputfile "\n boot obs median: $boot_obs_median boot exp median: $boot_exp_median \n";
			print $outputfile "\n boot obs mean: $boot_obs_mean boot exp mean: $boot_exp_mean \n";
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
			print $outputfile "bin $j mean_boot_obs $mean_obs mean_boot_exp $mean_exp diff_percentile_5 ".$stat_obs->percentile(5)." diff_percentile_95 ".$stat_obs->percentile(95).".\n";
		}
	
	
	
		close CSVFILE;
		print $outputfile "Number of iterations: $iteration\n";
		print $outputfile "- pvalue_epistasis  pvalue_environment\n";
		print $outputfile "median_stat ".($pval_epi/$iteration)." ".($pval_env/$iteration)."\n";
		print $outputfile "mean_stat ".($pval_epi_for_mean/$iteration)." ".($pval_env_for_mean/$iteration)."\n";
		$self->printFooter($outputfile);
		close $outputfile;	
		
		}
		else {
			print $outputfile "hist sum is 0";
			close $outputfile;
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
			## copyaste from prepare: create restricted hash	
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
			#	print "MAXBIN $maxbin\n";
			}
			
			## end f copypaste	
			my $count = scalar keys %obs_hash_restricted;
			print  COUNTER "$restriction ".$group_names[$group_number]." group $count "; 
			if ($count == 0){
				$group_number++;
				next;
			}
			my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]);
			my $outputfile;
			if ($fake){
				open $outputfile, ">>$file" or die "Cannot create $file";
			}
			else {
				open $outputfile, ">$file" or die "Cannot create $file";
			}
			
			
			#copypaste from all
			my %histhash;
			foreach my $site_node(keys %obs_hash_restricted){
				foreach my $bin (keys %{$obs_hash_restricted{$site_node}}){
					$histhash{$bin}[0] += $obs_hash_restricted{$site_node}{$bin}[0];
					$histhash{$bin}[1] += $obs_hash_restricted{$site_node}{$bin}[1];
				}
			}
			print $outputfile "bin\tobs\texp\n";
			my @sorted_bins = sort { $a <=> $b } keys %histhash;
			foreach my $bin (@sorted_bins){
				print $outputfile $bin."\t".$histhash{$bin}[0]."\t".$histhash{$bin}[1]."\n";
			}
			## end of copypaste
			
			
			my %flat_obs_hash;
			my %flat_exp_hash;
			#print  going to flat hash\n";
			foreach my $bin(1..$maxbin){
				foreach my $node (keys %obs_hash_restricted){
					$flat_obs_hash{$bin} += $obs_hash_restricted{$node}{$bin}[0];
					$flat_exp_hash{$bin} += $obs_hash_restricted{$node}{$bin}[1];
				}
			}
			#print  computing hist emdian \n"	;
			my $obs_median = hist_median_for_hash(\%flat_obs_hash, $step);
			my $exp_median = hist_median_for_hash(\%flat_exp_hash, $step);
			my $obs_mean = hist_mean_for_hash(\%flat_obs_hash, $step);
			my $exp_mean = hist_mean_for_hash(\%flat_exp_hash, $step);
			
			print $outputfile "\n observed median: $obs_median expected poisson median $exp_median observed mean: $obs_mean expected poisson mean $exp_mean\n";
			if($obs_mean eq "NaN" || $exp_mean eq "NaN") {
				print $outputfile " hist sum is 0";
				close $outputfile;
				$group_number++; # to skip complement for this group
				next;
			}
			#print going to read input file\n";		
				
			my $csvfile = File::Spec->catfile($outdir,temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number].".csv");
			open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
			my $iteration = 0; # counter for meaningful iterations
			my @group_boot_medians;
			my @group_boot_means;
			my @hist_obs;
			my @hist_exp;
			my @array_gbo_minus_gbe;
			my $itnumber = 0; # tracks iteration number, so that group and its complement are taken from the same iteration
			while(<CSVFILE>){
				my %boot_obs_hash;
				my %boot_exp_hash;
				my @splitter = split(/,/, $_);
				if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
					$itnumber++;
					next;
				}
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
					#push @{$array_gbo_minus_gbe[$bin10]}, $gbo_minus_gbe_10bin_array[$bin10];
					$array_gbo_minus_gbe[$bin10][$itnumber] = $gbo_minus_gbe_10bin_array[$bin10];
				}
				## end of copypaste
				
				
				## 5.02 redundant
				#for (my $i = 0; $i < scalar @splitter; $i++){
				#	my $bin = ($i/2)+1;
				#	$boot_obs_hash{$bin} = $splitter[$i];
				#	$i++;
				#	$boot_exp_hash{$bin} = $splitter[$i];
				#}
				
				my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash, $step);
				my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash, $step);
				my $boot_obs_mean = hist_mean_for_hash(\%boot_obs_hash, $step);
				my $boot_exp_mean = hist_mean_for_hash(\%boot_exp_hash, $step);
				
				$group_boot_medians[$itnumber][0] = $boot_obs_median;
				$group_boot_medians[$itnumber][1] = $boot_exp_median;
				$group_boot_means[$itnumber][0] = $boot_obs_mean;
				$group_boot_means[$itnumber][1] = $boot_exp_mean;			
				unless ($fake) {print $outputfile "\n boot obs median: $boot_obs_median boot exp median $boot_exp_median  boot obs mean: $boot_obs_mean boot exp mean $boot_exp_mean\n";}
				$itnumber++;
				$iteration++;
			}
			close CSVFILE;
			$self->printFooter($outputfile);	
			close $outputfile;
			
				print "Number of meaningful iterations for group ".$group_names[$group_number]." is $iteration (haven't looked at the complement yet)\n";	
				if ($iteration == 0) {
					$group_number++;
					next;
				}
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
			#	print "MAXBIN $maxbin\n";
			}
				
			## end of copypaste	
				
			my $count = scalar keys %complement_obs_hash_restricted;
			print  COUNTER " $restriction ".$group_names[$group_number]." complement $count\n"; 
			
			my $file = File::Spec->catfile($outdir,$prot."_gulpselector_vector_boot_median_test_".$restriction."_".$group_names[$group_number]);
			my $outputfile;
			if ($fake){
				open $outputfile, ">>$file" or die "Cannot create $file";
			}
			else {
				open $outputfile, ">$file" or die "Cannot create $file";
			}
			my %complement_flat_obs_hash;
			my %complement_flat_exp_hash;
			#print " going to flat hash\n";
			foreach my $bin(1..$maxbin){
				foreach my $node (keys %complement_obs_hash_restricted){
					$complement_flat_obs_hash{$bin} += $complement_obs_hash_restricted{$node}{$bin}[0];
					$complement_flat_exp_hash{$bin} += $complement_obs_hash_restricted{$node}{$bin}[1];
				}
			}
			#print " compung hist emdian \n"	;
			my $complement_obs_median = hist_median_for_hash(\%complement_flat_obs_hash, $step);
			my $complement_exp_median = hist_median_for_hash(\%complement_flat_exp_hash, $step);
			my $complement_obs_mean = hist_mean_for_hash(\%complement_flat_obs_hash, $step);
			my $complement_exp_mean = hist_mean_for_hash(\%complement_flat_exp_hash, $step);
			
			print $outputfile "\n observed median: $complement_obs_median expected median: $complement_exp_median observed mean: $complement_obs_mean expected mean: $complement_exp_mean\n";
			if($complement_obs_mean eq "NaN" || $complement_exp_mean eq "NaN") {
				print $outputfile " hist sum is 0";
				close $outputfile;
				next;
			}
			
			my $csvfile = File::Spec->catfile($outdir,temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$group_names[$group_number].".csv");
			open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
			my $iteration = 0;
			my @complement_boot_medians;
			my @complement_boot_means;
			my @hist_compl_obs;
			my @hist_compl_exp;
			my @array_cbo_minus_cbe;
			my $itnumber = 0;
			while(<CSVFILE>){
				my %boot_obs_hash;
				my %boot_exp_hash;
				my @splitter = split(/,/, $_);
				if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
					$itnumber++;
					next;
				}
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
					#push @{$array_cbo_minus_cbe[$bin10]}, $cbo_minus_cbe_10bin_array[$bin10];
					$array_cbo_minus_cbe[$bin10][$itnumber] = $cbo_minus_cbe_10bin_array[$bin10];
				}
				## end of copypaste
				
				
				## 5.02 redundant
				#for (my $i = 0; $i < scalar @splitter; $i++){
				#	my $bin = ($i/2)+1;
				#	$boot_obs_hash{$bin} = $splitter[$i];
				#	$i++;
				#	$boot_exp_hash{$bin} = $splitter[$i];
				#}
				
				my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash, $step);
				my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash, $step);
				my $boot_obs_mean = hist_mean_for_hash(\%boot_obs_hash, $step);
				my $boot_exp_mean = hist_mean_for_hash(\%boot_exp_hash, $step);			
				$complement_boot_medians[$itnumber][0] = $boot_obs_median;
				$complement_boot_medians[$itnumber][1] = $boot_exp_median;
				$complement_boot_means[$itnumber][0] = $boot_obs_mean;
				$complement_boot_means[$itnumber][1] = $boot_exp_mean;
				unless ($fake) {print $outputfile "\n boot obs median: $boot_obs_median boot exp median $boot_exp_median  boot obs mean: $boot_obs_mean boot exp mean $boot_exp_mean\n";}
				$itnumber++;
				$iteration++;
			}
			close CSVFILE;
		
			## 20.09.2016
			# delete iteration data from group, if there were no nodes in complement in this iteration, and vice versa		
			for (my $it = 0; $it < $itnumber; $it++){
				if (! defined $complement_boot_medians[$it] || ! defined $group_boot_medians[$it]){
					$complement_boot_medians[$it] = undef;
					$complement_boot_means[$it] = undef;
					for (my $bin10 = 0; $bin10 < scalar @array_cbo_minus_cbe; $bin10++){
						$array_cbo_minus_cbe[$bin10][$it] = undef;
					}
					$group_boot_medians[$it] = undef;
					$group_boot_means[$it] = undef;
					for (my $bin10 = 0; $bin10 < scalar @array_gbo_minus_gbe; $bin10++){
						$array_gbo_minus_gbe[$bin10][$it] = undef;
					}
				} 
			}
			# now group and complement contain undefs at the same indices
			@complement_boot_medians = grep defined, @complement_boot_medians;
			@complement_boot_means = grep defined, @complement_boot_means;
			@group_boot_medians = grep defined, @group_boot_medians;
			@group_boot_means = grep defined, @group_boot_means;
			for (my $bin10 = 0; $bin10 < scalar @array_cbo_minus_cbe; $bin10++){
					@{$array_cbo_minus_cbe[$bin10]} = grep defined, @{$array_cbo_minus_cbe[$bin10]};
			}
			for (my $bin10 = 0; $bin10 < scalar @array_gbo_minus_gbe; $bin10++){
					@{$array_gbo_minus_gbe[$bin10]} = grep defined, @{$array_gbo_minus_gbe[$bin10]};
			}
			##
			print "Number of meaningful iterations (used for pvalue estimation) for group ".$group_names[$group_number]." is ".scalar @complement_boot_medians." or ".scalar @group_boot_means." first one is used for division, second - for iterating through simulations\n";
			
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
			
			if (scalar @complement_boot_medians != scalar @group_boot_medians){
				print "Error! complement_boot_medians size is not equal to group_boot_medians size: ".scalar @complement_boot_medians." != ".scalar @group_boot_medians."\n";
			}
			my $updated_iteration_number = scalar @complement_boot_medians;
			
			unless($fake){
			print $outputfile "bin mean_group_boot_obs mean_group_boot_exp ";
			print $outputfile "mean_compl_boot_obs mean_compl_boot_exp ";
			print $outputfile " diff_gbo-gbe-cbo+cbe_percentile_5 diff_gbo-gbe-cbo+cbe_percentile_95";
			print $outputfile " diff_gbo-gbe_percentile_5 diff_gbo-gbe_percentile_95";
			print $outputfile " diff_cbo-cbe_percentile_5 diff_cbo-cbe_percentile_95\n";
			my $maxbin  = max(scalar @hist_obs, scalar @hist_compl_obs);
			$maxbin = max ($maxbin, scalar @hist_exp);
			$maxbin = max ($maxbin, scalar @hist_compl_exp);
			
			
			for (my $j = 0; $j < $maxbin; $j++){ #foreach bin
				my $mean_group_obs = $hist_obs[$j]/$updated_iteration_number; # 26.09.2016 - $updated_iteration_number instead of $iteration
				my $mean_group_exp = $hist_exp[$j]/$updated_iteration_number;
				my $mean_compl_obs = $hist_compl_obs[$j]/$updated_iteration_number;
				my $mean_compl_exp = $hist_compl_exp[$j]/$updated_iteration_number;
				my $stat_gbo_minus_gbe = Statistics::Descriptive::Full->new();
				$stat_gbo_minus_gbe->add_data(\@{$array_gbo_minus_gbe[$j]});
				my $stat_cbo_minus_cbe = Statistics::Descriptive::Full->new();
				$stat_cbo_minus_cbe->add_data(\@{$array_cbo_minus_cbe[$j]});	
				my $stat_diffdiff = Statistics::Descriptive::Full->new();
				$stat_diffdiff->add_data(\@{$array_diffdiff[$j]});			
				print $outputfile "$j $mean_group_obs $mean_group_exp ";
				print $outputfile "$mean_compl_obs $mean_compl_exp ";
				print $outputfile $stat_diffdiff->percentile(5)." ".$stat_diffdiff->percentile(95);
				print $outputfile $stat_gbo_minus_gbe->percentile(5)." ".$stat_gbo_minus_gbe->percentile(95);
				print $outputfile $stat_cbo_minus_cbe->percentile(5)." ".$stat_cbo_minus_cbe->percentile(95)."\n";
			}
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
			print $outputfile "Number of iterations: ".$updated_iteration_number."\n"; # and not $iteration. changed at 26.09.2016 
			print $outputfile "- pvalue_epistasis_enrichment pvalue_environment_enrichment pvalue_epistasis pvalue_environment\n";
			print $outputfile "median_stat ".($pval_epi_enrichment/$updated_iteration_number)." ".($pval_env_enrichment/$updated_iteration_number)." ".($pval_epi/$updated_iteration_number)." ".($pval_env/$updated_iteration_number)."\n";
			print $outputfile "mean_stat ".($pval_epi_enrichment_for_mean/$updated_iteration_number)." ".($pval_env_enrichment_for_mean/$updated_iteration_number)." ".($pval_epi_for_mean/$updated_iteration_number)." ".($pval_env_for_mean/$updated_iteration_number)."\n";
			$self->printFooter($outputfile);
			close $outputfile;	
			
			
		}
		
	}
	close COUNTER;
}

# 13.09.2016 pvalues for every ancestor node in the tree	
sub count_single_site_pvalues{	
	my $self = $_[0];
	my $prot = $self -> {static_protein};
	my @restriction_levels = @{$_[1]};
	my $dir = $self -> {static_output_base};
	my $outdir = $self -> {static_output_subfolder};
	
	my $realdata =  $self -> {realdata};
	my $maxbin = $realdata->{"maxbin"}; 
	my $step = $realdata->{"step"}; #bin size
	unless (defined $step) {die "Oh no, bin size in realdata is not defined. Won't proceed with simulations.\n";}
	#print "before cycle\n";
	
	my $obs_hash = get_obshash($realdata, List::Util::min(@restriction_levels)); # if min $restriction is less than restriction in realdata, it will die
	my $subtree_info = $realdata->{"subtree_info"};

	#my $restriction = List::Util::min(@restriction_levels); # todo: print maxdepth instead of restriction level in filenames (or somewhere inside the output file), so that we analyse each site_node only once
	for my $restriction(@restriction_levels){
		print "level $restriction\n";
		my %obs_hash_restricted;
		my %norm_restricted;
		
		## create restricted hash	
		foreach my $site_node(keys %{$obs_hash}){
			my ($site, $node_name) = split(/_/, $site_node);
			my $maxdepth = $subtree_info->{$node_name}->{$site}->{"maxdepth"};
			if ($maxdepth > $restriction){ 
			foreach my $bin(keys %{$obs_hash->{$site_node}}){
					$maxbin = max($bin, $maxbin);
					$norm_restricted{$site_node} += $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[0] = $obs_hash->{$site_node}->{$bin}->[0];
					$obs_hash_restricted{$site_node}{$bin}[1] = $obs_hash->{$site_node}->{$bin}->[1];
				}
		
			}
		}
		##
		my $count = scalar keys %obs_hash_restricted;
	
		my $file = File::Spec->catfile($outdir, $prot."_gulpselector_vector_boot_median_test_".$restriction."_single_sites");
		open my $outputfile, ">$file" or die "Cannot create $file";


		foreach my $site_node (keys %obs_hash_restricted){
			my %flat_obs_hash;
			my %flat_exp_hash;
			foreach my $bin(1..$maxbin){
				$flat_obs_hash{$bin} += $obs_hash_restricted{$site_node}{$bin}[0]; 
				$flat_exp_hash{$bin} += $obs_hash_restricted{$site_node}{$bin}[1]; 
			}
			my $obs_median = hist_median_for_hash(\%flat_obs_hash, $step);
			my $exp_median = hist_median_for_hash(\%flat_exp_hash, $step);
			my $obs_mean = hist_mean_for_hash(\%flat_obs_hash, $step); # 18.03 - added the same statistics based on histogram mean (instead of median)
			my $exp_mean = hist_mean_for_hash(\%flat_exp_hash, $step);
			

			print $outputfile "site_node\tbin\tobs\texp\n";
				my @sorted_bins = sort { $a <=> $b } keys $obs_hash_restricted{$site_node};
				foreach my $bin (@sorted_bins){
					if (defined $obs_hash_restricted{$site_node}{$bin}[0] && defined $obs_hash_restricted{$site_node}{$bin}[1]){
						print $outputfile $site_node."\t".$bin."\t".$obs_hash_restricted{$site_node}{$bin}[0]."\t".$obs_hash_restricted{$site_node}{$bin}[1]."\n";
					}
				}

			
			
			print $outputfile "\n site_node: $site_node\n";
			print $outputfile "\n observed median: $obs_median\n";
			print $outputfile "\n poisson expected median: $exp_median\n";
			print $outputfile "\n observed mean: $obs_mean\n";
			print $outputfile "\n poisson expected mean: $exp_mean\n";
			
			my $pval_epi;
			my $pval_env;
			my $pval_epi_for_mean;
			my $pval_env_for_mean;
			
		
		if ($obs_mean ne "NaN"){
		
		my $csvfile = File::Spec->catfile($outdir, temp_tag(),$prot."_gulpselector_vector_".$restriction."_".$site_node.".csv");
		open CSVFILE, "<$csvfile" or die "Cannot open $csvfile";
		my $iteration = 0;
		my @hist_obs;
		my @hist_exp;
		my @array_obs_minus_exp;
		while(<CSVFILE>){
			my %boot_obs_hash;
			my %boot_exp_hash;
			my @splitter = split(/,/, $_);
			if ($splitter[0] eq "NA"){ # if no appropriate nodes were produced in this iteration, it is skipped
				next;
			}
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
			
			my $boot_obs_median = hist_median_for_hash(\%boot_obs_hash, $step);
			my $boot_exp_median = hist_median_for_hash(\%boot_exp_hash, $step);
			my $boot_obs_mean = hist_mean_for_hash(\%boot_obs_hash, $step);
			my $boot_exp_mean = hist_mean_for_hash(\%boot_exp_hash, $step);
			print $outputfile "\n boot obs median: $boot_obs_median boot exp median: $boot_exp_median \n";
			print $outputfile "\n boot obs mean: $boot_obs_mean boot exp mean: $boot_exp_mean \n";
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
			print $outputfile "bin $j mean_boot_obs $mean_obs mean_boot_exp $mean_exp diff_percentile_5 ".$stat_obs->percentile(5)." diff_percentile_95 ".$stat_obs->percentile(95).".\n";
		}
	
		close CSVFILE;
		
		#print FILE "- pvalue_epistasis  pvalue_environment\n";
		#print FILE "median_stat ".($pval_epi/$iteration)." ".($pval_env/$iteration)."\n";
		#print FILE "mean_stat ".($pval_epi_for_mean/$iteration)." ".($pval_env_for_mean/$iteration)."\n";
		my ($psite, $pnode_name) = split(/_/, $site_node);
		my $pmaxdepth = $subtree_info->{$pnode_name}->{$psite}->{"maxdepth"};
		my $pmutcount = sum(values %flat_obs_hash);
		print $outputfile "Number of iterations: $iteration\n";
		print $outputfile "#\tsite_node\tmutations\tmaxlength\tpvalue_epistasis(median)\tpvalue_epistasis(mean)\tpvalue_environment(median)\tpvalue_environment(mean)\n";
		print $outputfile ">\t".$site_node."\t".$pmutcount."\t".$pmaxdepth."\t".($pval_epi/$iteration)."\t".($pval_epi_for_mean/$iteration)."\t".($pval_env/$iteration)."\t".($pval_env_for_mean/$iteration)."\n";
		}
		else {
			print $outputfile "hist sum is 0";	
		}
		}
			$self -> printFooter($outputfile);
		close $outputfile;	
	}	

}


sub no_check{
	return 1;
}


# prints protein_for_LRT files
sub print_data_for_LRT {
	my $self = shift;
	# my $dir = File::Spec->catdir(getcwd(),"likelihood", $self->{static_state}); # before august 2016 refactoring 
	my $dir = File::Spec -> catfile($self->{static_output_base}, "likelihood");
	make_path($dir);
	my $filename = File::Spec->catfile($dir, ($self->{static_protein})."_for_LRT.csv");
	open my $file, ">$filename" or die "Cannot create $filename";
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my @group = (1..$self->mylength());
	
	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ($root);
	$self->visitor_coat ($root, \@array,\&lrt_visitor,\&no_check,\@args,0);
	print $file "site,ancestor_node,t_branch_start,t_branch_end,event_indicator\n";
	foreach my $ind (@group){
		foreach my $ancnode(@{$self->{static_nodes_with_sub}{$ind}}){
			if(ref($ancnode) eq "REF"){
				$ancnode = ${$ancnode};
			}
			my $ancnodename = $ancnode->get_name();
			foreach my $node ( keys %{$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}}){
				print $file $ind.",".$ancnodename.",".$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[1].",".
				$self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[2].",";
				my $event = 0;
				if($self->{static_subtree_info}{$ancnodename}{$ind}{"lrt"}{$node}[0]) {$event = 1};
				print $file "$event\n";
			}
		}
	}
	$self->printFooter($file);	
	close $file;
}




## Since 5.02
sub depth_groups_entrenchment_optimized_selection_alldepths {
	my $self = shift;
	my $step = $_[0];
	my $restriction = $_[1]; # 10.10 restriction returns
	my $ancestral_nodes = $_[2];
	my $overwrite = $_[3];
	my $tag = $_[4];
	my $verbose = $_[5];
	my $gr = $_[6];
	my $dir = File::Spec->catdir($self->{static_output_base}, $self->{static_protein});
	make_path($dir);
	my $filename = File::Spec->catfile($dir, $self->{static_protein}."_for_enrichment_".$tag);
	if ($verbose){print "Going to print in file $filename \n";}
	open my $file, ">>$filename" or die "Cannot create $filename\n";
	my $root = $self->{static_tree}-> get_root;
	my @array;
	my %hist;
	print $file ">new iteration\n";
	
	my @group;
	if ($gr){ 
		@group = @{$gr};
	}
	else {
		my $length = $self->mylength();
		@group = (1..$length);
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
		foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node;
			if(ref($nod) eq "REF"){
				$node = ${$nod};
			}
			else {$node = $nod};
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
			if ($self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > $restriction){ #10.10 restriction returns
				
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
					print $file "site $ind node ".$node->get_name()." maxdepth ".$self->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}."\n";
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

	 print $file "$bin,".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";
	 $check_total_obs += $hist{$site_node}{$bin}[0];
	 $check_total_exp += $hist{$site_node}{$bin}[1];
				}

				if ($total_length == $check_local_lengths_sum){
				#print "local lengths sumtest ok: $total_length $check_local_lengths_sum\n";
				}
				else {
				print "local length sumtest failed! total $total_length, local_sum $check_local_lengths_sum\n";
				}
				if ($check_total_obs-$check_total_exp < 0.001 && -$check_total_obs+$check_total_exp < 0.001 ){
				#print "obsexp sumtest ok\n";
				}
				else {
				print "obsexp sumtest failed! total obs $check_total_obs, total exp $check_total_exp total_muts $total_muts site $ind node ".$node->get_name()."\n";
				}
				}

			}
			
		}
	}	
	close $file;
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
		@group = (1..$self->mylength());
	}
	

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root);
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	my $debugnum = scalar keys %{$self ->{static_nodes_with_sub}};
	print "News from depth..2: static_nodes_with_sub contains $debugnum keys\n";
	my $debugnum = scalar keys %{$self ->{static_subtree_info}};
	print 	"News from depth..2: static_subtree_info contains $debugnum keys (nodes)\n";
	foreach my $ind (@group){
		foreach my $nod(@{$self ->{static_nodes_with_sub}{$ind}}){
			print "nod is ".$nod." ref(nod) is ".ref($nod)."\n";
			my $node;
			if(ref($nod) eq "REF"){
				$node = ${$nod};
			}
			else {$node = $nod;}
			my $site_node = $ind."_".$node->get_name();
			print "site_node $site_node \n";
			my $total_muts;
			my $total_length;
	#print "depth ".$static_depth_hash{$ind}{$node->get_name()}."\n";
			print " maxdepth is ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}. " and restriction is $restriction\n";
			if ($self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"} > $restriction){
				print "will try it\n";
				my %subtract_hash;
				
				if ($self ->{static_subtract_tallest}){
					print "just checking ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth_node"}."\n";
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
				print $node->get_name()." ".$ind." TOTALS: $total_muts, $total_length, maxdepth ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}."\n";
				#if ($total_length > 0 && $total_muts/$total_length < 0.005){
				if ($total_length > 0 && $total_muts > 0){ # 21.10 added total_muts > 0
				print "total muts $total_muts \n";
			#	print "site $ind node ".$node->get_name()." maxdepth ".$self ->{static_subtree_info}{$node->get_name()}{$ind}{"maxdepth"}."\n"; # commented out 16.09
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
	# print "$bin,".$hist{$site_node}{$bin}[0].",".$hist{$site_node}{$bin}[1]."\n";  # commented out 16.09
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
						my $str =  $ind."_".$node->get_name();
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
	$self->visitor_coat ($root, \@array,\&synresearch_visitor,\&no_check,\@args,0);
	foreach my $ind (@group){
			foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node = ${$nod};
			if (!$node->is_terminal){
				my $anc = $self->{static_subs_on_node}{$node->get_name()}{$ind}->{"Substitution::derived_allele"};
			#	print "anc $anc\n";
				my $synmuts = compare::get_synmuts($anc);
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
		@group = (1..$self->mylength());
	}
	 my $name = $_[4];

	my %closest_ancestors;
	$root->set_generic("-closest_ancestors" => \%closest_ancestors);
	my @args = ( $step, $root);
	$self->visitor_coat ($root, \@array,\&entrenchment_visitor,\&no_check,\@args,0);
	my %sitecounts;
	my %mutcounts;
	foreach my $ind (@group){
		foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node;
			if(ref($nod) eq "REF"){
				$node = ${$nod};
			}
			else {$node = $nod;}
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


#takes a hash of probabilities for 0,1,2...
sub hist_median_for_hash{
	my @hist = hist_to_array($_[0]);
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	if (! defined $step) {$step  = 1;}
	return hist_median(\@hist, $step);
}

#takes a hash of probabilities for 0,1,2...
sub hist_mean_for_hash {
	my @hist = hist_to_array($_[0]);
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	if (! defined $step) {$step  = 1;}
	return hist_mean(\@hist, $step);
}

#takes an array of probabilities for 0,1,2...
sub hist_median{
	my @hist = @{$_[0]};
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	if (! defined $step) {$step  = 1;}
	my $summ = sum (@hist);
	my $head = 0;
	my $interval = 0;
	my $median = 0;
	
	while ($head < $summ/2){
		$head += $hist[$interval];
		$median = $interval*$step; 
		$interval++;
	}
	
	if ($head == $summ/2){
	#	$median += 0.5*$step;
		my $leftmedian = $median;
		my $rightmedian;
		my $newhead = $head;
		while ($newhead == $head){
			$newhead += $hist[$interval];
			$rightmedian = $interval*$step; 
			$interval++;
		}
		$median = ($leftmedian+$rightmedian)/2;
	}
#print_hist(\@hist);
	return $median;
}

sub hist_mean {
	my @hist = @{$_[0]};
	my $step = $_[1]; # 21.09.2016 for bin size = 0.5
	if (! defined $step) {$step  = 1;}
	my $summ = sum (@hist);
	my $integer;
	for(my $i = 0; $i <scalar @hist; $i++){
		$integer += $i*$hist[$i]*$step;
	}
	if ($summ > 0){
		return $integer/$summ;
	}
	else {
		return "NaN";
	}
}



# takes hash of probabilities for 0,1,2... and returns an array of ordered values
sub hist_to_array {
	my %prehist =  %{$_[0]};
	my @hist;
	my @sorted_keys = sort {$a <=> $b} keys %prehist;
	# for (my $i = 1; $i <= $sorted_keys[-1]; $i++){ # changed at 21.09.2016
	for (my $i = 0; $i <= $sorted_keys[-1]; $i++){
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

# one hist for one ancestor aa in one site
# the chosen one (used by r scripts for drawing plots)
# circles, not rings! ()
# mutations at the ends of branches !

sub egor_smart_site_entrenchment {
	my $self = shift;
	my $verbose = shift;
	my $step = 1;
	my $root = $self->{static_tree}-> get_root;
	my $file = File::Spec->catfile($self->{static_output_base}, "egor_smart_".$self->{static_protein}.".csv");
	open my $plotcsv, ">$file" or die "Cannot create $file \n";
	my @array;
	print $plotcsv "radius,site,node,density,cummuts,cumlength\n";
	my $hash_ready;
	if (exists $self->{static_ring_hash}){
		warn "Static_ring_hash is ready, egor_smart_site_entrenchment won't change it\n";
		$hash_ready = 1;
	}

	for (my $ind = 1; $ind < $self->mylength(); $ind++){
		if ($verbose) {print "$ind\n"};
		foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node = ${$nod};
			if (!$node->is_terminal){
			if ($verbose) {print $node->get_name()."\n"};
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			unless ($hash_ready) {$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);} #visitor_coat cannot be used inside a cycle; we still do not want to mess the hash up, so we check for its existance before this cycle
			
			my $cumulative_muts;
			my $cumulative_length;
			my @sorted_keys = sort {$a <=> $b} (keys %{$self->{static_ring_hash}{$ind}{$node->get_name()}});
			foreach my $bin (@sorted_keys){
				my $muts_in_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[1];
				$cumulative_muts += $muts_in_bin;
				$cumulative_length += $length_of_bin;
				if ($verbose) {print "bin $bin observed ".$self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0]." totmut $cumulative_muts muts in bin $muts_in_bin totlen $cumulative_length bin len $length_of_bin\n"};
				if ($cumulative_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0]/$self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[1]; #density

					$hist{$bin} = $cumulative_muts/$cumulative_length; #density
				}
				else {
					$hist{$bin} = 0;
				}
					if ($muts_in_bin > 0){	
					print $plotcsv "$bin,$ind,".$node->get_name().",".$hist{$bin}.",".$cumulative_muts.",".$cumulative_length."\n";
				}
				
			}
			

		}
		}
		
	}

	$self->printFooter($plotcsv);
	close $plotcsv;
	
}

# last egor plots sub. NOT the chosen one (egor_h1.csv are printed by some other method (discovered by comparison))
# honest rings, can be used for mann-kendall analysis
# ! mutations at the ends of branches !

sub egor_diff_rings_site_entrenchment {
	my $self = shift;
	print ($self->{static_protein}."\n");
	my $step = 1;
	my $root = $self->{static_tree}-> get_root;
	my $file = File::Spec->catfile($self->{static_output_base}, "egor_diff_rings_".$self->{static_protein}.".csv");
	open my $plotcsv, ">$file" or die "Cannot create $file \n";
	my @array;
	print $plotcsv "radius,site,node,density,cum_muts,cum_length\n";
	
	my $hash_ready;
	if (exists $self->{static_ring_hash}){
		$hash_ready = 1;
		warn "Static_ring_hash is ready, egor_diff_rings_site_entrenchment won't change it\n";
	}

	for (my $ind = 1; $ind < $self->mylength(); $ind++){
		#print "$ind\n";
		foreach my $nod(@{$self->{static_nodes_with_sub}{$ind}}){
			my $node = ${$nod};
			if (!$node->is_terminal){
			#print $node->get_name()."\n";
			my %hist;
			#my $ancestor = ${$static_subs_on_node{$node->get_name()}}{$ind}->{"Substitution::ancestral_allele"};
			my @args = ($ind, $step, $node);
			unless ($hash_ready) {$self->my_visit_depth_first ($node, \@array,\&update_ring,\&has_no_mutation,\@args,0);} #visitor_coat cannot be used inside a cycle; we still do not want to mess the hash up, so we check for its existance
			my $cumulative_muts;
			my $cumulative_length;
			my @sorted_keys = sort {$a <=> $b} (keys %{$self->{static_ring_hash}{$ind}{$node->get_name()}});
			foreach my $bin (@sorted_keys){
				#print "bin $bin observed ".$self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0]." totmut $cumulative_muts totlen $cumulative_length\n";
				my $muts_in_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[0];
				my $length_of_bin = $self->{static_ring_hash}{$ind}{$node->get_name()}{$bin}[1];
				$cumulative_muts += $muts_in_bin;
				$cumulative_length += $length_of_bin;
				if ($cumulative_length > 0){ #there are some internal nodes with 0-length terminal daughter branches
					#	$hist{$bin}{$ancestor} += $static_ring_hash{$ind}{$node->get_name()}{$bin}[0]/$static_ring_hash{$ind}{$node->get_name()}{$bin}[1]; #density
					$hist{$bin} = $cumulative_muts/$cumulative_length; #density
				}
				else {
					$hist{$bin} = 0;
				}
				#if ($muts_in_bin > 0 || $bin == $sorted_keys[0] || $bin == $sorted_keys[-1]){	
					if ($muts_in_bin > 0){	
					print $plotcsv "$bin,$ind,".$node->get_name().",".$hist{$bin}.",".$cumulative_muts.",".$cumulative_length."\n";
					$cumulative_muts = 0;
					$cumulative_length = 0;
				}
				
			}
			

		}
		}
	}
	
	close $plotcsv;
	
	
}


# decorator for my_visit_depth_first, checks for existance of corresponding hash and prevents unintentional changes in it (or deletes it and overwrites)
# must not be used in loop context!
## !!! does not work as expected, should not be used for preventing changes in pre-existing data. 
sub visitor_coat {
		my $self = shift;
		my $node = $_[0];
		my @array = @{$_[1]};
		my $action_callback = $_[2];
		my $check_callback = $_[3];
		my $callback_args = $_[4];
		my $depth = $_[5];
		my $overwrite = $_[6];
		
		my $visitor_name = sub_name($action_callback);
		
		if ($visitor_name eq "update_ring"){
			if (exists $self->{static_ring_hash}){
				if ($overwrite){
					delete $self->{static_ring_hash};
				}
				else {return;}	
			}	
		}
		elsif ($visitor_name eq "entrenchment_visitor" || $visitor_name eq "lrt_visitor"){
			if (exists $self->{static_subtree_info}){
				if ($overwrite){
					delete $self->{static_subtree_info};
				}
			# commented out at 27.09.2016	
			#		foreach my $node(keys $self->{static_subtree_info}){
			#				foreach my $site(keys $self->{static_subtree_info}{$node}){
			#					if ($visitor_name eq "entrenchment_visitor"){
			#						if (exists $self->{static_subtree_info}{$node}{$site}{"hash"}){
			#							if ($overwrite){
		#									delete $self->{static_subtree_info}{$node}{$site}{"hash"};
	#										delete $self->{static_subtree_info}{$node}{$site}{"maxdepth"};
	#										delete $self->{static_subtree_info}{$node}{$site}{"maxdepth_node"};
	#									}
	#									else {return;}
	#								}
	#									
	#							}
	#							elsif ($visitor_name eq "lrt_visitor"){
	#								if (exists $self->{static_subtree_info}{$node}{$site}{"lrt"}){
	#									if ($overwrite){
	#										delete $self->{static_subtree_info}{$node}{$site}{"lrt"};
	#									}
	#									else {return;}
	#								}
	#							}
	#						}
	#				}
			}
		}
		else {print "Warning: visitor_coat cannot perform any check for $visitor_name subroutine. Launching my_visit_depth_first without checking for previous launches\n";}
		$self->my_visit_depth_first($node, \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
}

 sub my_visit_depth_first {
		my $self = shift;
		my $node = $_[0];
		my @array = @{$_[1]};
		my $action_callback = $_[2];
		my $check_callback = $_[3];
		my $callback_args = $_[4];
		my $depth = $_[5];
		push @array, $node;
		my $len = $node -> get_branch_length;
		$depth += $len;
		&$action_callback($self, $node, $callback_args, $depth);
		if (! $node->is_terminal && &$check_callback($self, $node, $callback_args)){
			my $i = 0;
			while($node->get_child($i)){
				@array = $self -> my_visit_depth_first($node->get_child($i), \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
				$i++;
				#@array = my_visit_depth_first($node->get_first_daughter, \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
				#@array = my_visit_depth_first($node->get_last_daughter, \@array, \&$action_callback, \&$check_callback, $callback_args, $depth);
			}
		}
		
		$node = pop @array;
#print $node->get_name()."\t".$depth."\n";
		$depth -= $len;
		return @array;
    }
    

    
    # old version, does not account for  mutations of the other type.
    # Used in update ring (for making plots) (that's ok, because plotting subroutines use newer has_no_mutation method to account for background mutations )
 	sub has_no_same_type_mutation { #has_no_same_type_mutation
 		my $self = $_[0];
 		my $node = $_[1];
 		my $site_index = $_[2]->[0];
 		my $starting_node = $_[2]->[2];
 		
 		if ($node eq $starting_node){
 			return 1;
 		}
 		if (${ $self->{static_subs_on_node}{$node->get_name()}}{$site_index}){
 			return 0;
 		}
 		else {
 			return 1;
 		}
 	}  
 	
 	
 	# also accounts for mutations of the other type (synonimous for nsyn and non-synonimous for syn)
 	sub has_no_mutation{
 		my $self = $_[0];
 		my $node = $_[1];
 		my $site_index = $_[2]->[0];
 		my $starting_node = $_[2]->[2];
 		
 		if ($node eq $starting_node){
 			return 1;
 		}
 		if (${$self->{static_subs_on_node}{$node->get_name()}}{$site_index}){
 			return 0;
 		}
 		
 		if ($self->{static_state} eq "nsyn"){
 			if (compare::is_neighbour_changing(${$self->{static_background_subs_on_node}{$node->get_name()}}{$site_index}, 1) == 1){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 		else {
 			if (${$self->{static_background_subs_on_node}{$node->get_name()}}{$site_index}){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 	}  
 	
 	
 	sub has_no_background_mutation {
 		my $self = $_[0];
 		my $node = $_[1];
 		my $site_index = $_[2];
 		if ($self->{static_state} eq "nsyn"){
 			if (compare::is_neighbour_changing(${$self->{static_background_subs_on_node}{$node->get_name()}}{$site_index}, 1) == 1){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 		else {
 			if (${$self->{static_background_subs_on_node}{$node->get_name()}}{$site_index}){
 				return 0;
 			}
 			else {
 				return 1;
 			}
 		}
 	}
 	

 	

 	
 	sub max ($$) { $_[$_[0] < $_[1]] }
 	
 	sub update_ring {
 		my $self = $_[0];
 		my $node = $_[1];
		my $site_index = $_[2]->[0];
		my $step = $_[2]->[1];
		my $starting_node = $_[2]->[2];
		my $depth = $_[3] - $starting_node->get_branch_length ;
		if ($node eq $starting_node){
			#print " \n equality: ".$starting_node ->get_name."\t".$node ->get_name."\n";
			return;
		}
		if ($starting_node -> is_terminal){
			#print " \n terminal: ".$starting_node ->get_name."\n";
			return;
		}
 		#print "depth $depth step $step bin ".(bin($depth,$step))."\n";
 		if (!($self->has_no_same_type_mutation($_[1], \@{$_[2]}))){ #has_no_same_type_mutation
 			$self->{static_ring_hash}{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[0] += 1;
 			#$self->{static_ring_hash}{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[1] -= ($node->get_branch_length)/2 # 19.09.2016 mutations happen in the middle of a branch; todo
 			#print "addded to 0\n";
 		}
 		#my $newdepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()} - ($node->get_branch_length)/2; # 19.09.2016 mutations happen in the middle of a branch; todo
 		$self->{static_ring_hash}{$site_index}{$starting_node->get_name()}{bin($depth,$step)}[1] += $node->get_branch_length;
 		#print "addded to 1\n";
 	}
 	
 	sub set_distance_matrix {
 		my $self = shift;
 		my $prot = $self->{static_protein};
 		my $file = File::Spec->catfile($self->{static_input_base}, $prot."_distance_matrix.csv");
 		if (! (-e $file)){
 				print "Preparing distance matrix..\n";
 				my $logs = capture ('Rscript Distances.R --treefile '. $self->{static_treefile}.' --output '.$file);
 				print $logs."\n";
 		}
	 	open CSV, "<$file" or die "Cannot open file $file\n";
	 	my $header = <CSV>;
	 	$header =~ s/[\s\n\r\t]+$//s;
	 	my @nodelables = split(',', $header);
	 	while(<CSV>){
				$_ =~ s/[\s\n\r\t]+$//s;
	 			my @dists = split(',', $_);
	 			my $node = $dists[0];
	 			for (my $i = 1; $i < scalar @dists; $i++){
	 				$self->{static_distance_hash}{$node}{$nodelables[$i]} = $dists[$i];
	 			}
	 	}
	 	close CSV;
 		

 	}
 	
 	
 	
 	
   # track_tallest is needed for finding longest path in the subtree and subtracting its length
   # Added at 08.10 for testing whether this will improve correspondence between simulation_observed and simulation_expected.
   # stopped using it at 15.09.2016    
   # 15.09.2016 version: halves of branches with foreground mutations are trimmed 
   # 19.09.2016: corrected
   	sub entrenchment_visitor {
 		my $self = shift;
 		my $node = $_[0];
 		my $step = $_[1]->[0];
 		my $subtract_tallest = $self->{static_subtract_tallest};
 		my $no_neighbour_changing = $self->{static_no_neighbour_changing};
 		my $no_leaves = $self->{static_no_leaves};
		if (!$node->is_root){
		my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")}; # closest_ancestors: ancestor mutation node for this node, key is a site number
		
		if (%closest_ancestors){
			foreach my $site_index(keys %closest_ancestors){ 
				my $anc_node = $closest_ancestors{$site_index};
				my $depth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()};
				$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($depth,$step)}[1] += $node->get_branch_length;
			#	print "anc ".$anc_node->get_name()." site ".$site_index." node ".$node->get_name()." depth $depth bin ".bin($depth,$step)." branchlength ".$node->get_branch_length."\n";
				my $current_maxdepth = $self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth"};
				if ($current_maxdepth < $depth){
						$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth"} = $depth;
						if ($subtract_tallest){
							$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth_node"} = $node->get_name();
						}
				}					
			}
			
			my @ancestors = keys %closest_ancestors;	
			foreach my $site_index(@ancestors){
				if (!($self->has_no_background_mutation($node, $site_index))){
					delete $closest_ancestors{$site_index};
				}
			}	
		}
		
		foreach my $site_index(keys %{$self->{static_subs_on_node}{$node->get_name()}}){
			
			if ($closest_ancestors{$site_index}){
				my $anc_node = $closest_ancestors{$site_index};
				my $halfdepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()} - ($node->get_branch_length)/2; #19.09.2016 
				my $fulldepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()}; #19.09.2016 
			#	print " ancestor ".$anc_node->get_name(). " node ".$node->get_name()." depth $depth\n";
			#	push $static_subtree_info{$anc_node->get_name()}{$site_index}{"nodes"}, \$node;
				if (!$no_neighbour_changing || ($no_neighbour_changing && ! compare::is_neighbour_changing($self->{static_subs_on_node}{$node->get_name()}{$site_index}, 1))){
					if (!$no_leaves || ($no_leaves && !($node->is_terminal()))){
						$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($halfdepth,$step)}[0] += 1; #19.09.2016 
					}
				}
				$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($halfdepth,$step)}[1] += ($node->get_branch_length)/2; # #19.09.2016  15.09.2016 version: halves of branches with foreground mutations are trimmed (the only thing I changed here) 
				$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($fulldepth,$step)}[1] -= $node->get_branch_length; #19.09.2016 we added this length before, but shouldn't have done it
			#	print "mutation! anc ".$anc_node->get_name()." site ".$site_index." node ".$node->get_name()." depth $depth bin ".bin($depth,$step)." branchlength ".$node->get_branch_length."\n";
				
			}
			$closest_ancestors{$site_index} = $node;
			
		}
		
		$node->set_generic("-closest_ancestors" => \%closest_ancestors);
		}
#!
 	}
  
  # entrenchment_visitorversion for for synmut_types: collect substitutions for each ancestor mutation (not just count them; do not track distance)
    	sub synresearch_visitor {
 		my $self = shift;
 		my $node = $_[0];
 		my $step = $_[1]->[0];
 		my $subtract_tallest = $self->{static_subtract_tallest};
 		my $no_neighbour_changing = $self->{static_no_neighbour_changing};
 		my $no_leaves = $self->{static_no_leaves};
		if (!$node->is_root){
		my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")}; # closest_ancestors: ancestor mutation node for this node, key is a site number
		
		if (%closest_ancestors){
			foreach my $site_index(keys %closest_ancestors){ 
				my $anc_node = $closest_ancestors{$site_index};
				my $depth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()};
				$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"hash"}{bin($depth,$step)}[1] += $node->get_branch_length;
			#	print "anc ".$anc_node->get_name()." site ".$site_index." node ".$node->get_name()." depth $depth bin ".bin($depth,$step)." branchlength ".$node->get_branch_length."\n";
				my $current_maxdepth = $self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth"};
				if ($current_maxdepth < $depth){
						$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth"} = $depth;
						if ($subtract_tallest){
							$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"maxdepth_node"} = $node->get_name();
						}
				}					
			}
			
			my @ancestors = keys %closest_ancestors;	
			foreach my $site_index(@ancestors){
				if (!($self->has_no_background_mutation($node, $site_index))){
					delete $closest_ancestors{$site_index};
				}
			}	
		}
		
		foreach my $site_index(keys %{$self->{static_subs_on_node}{$node->get_name()}}){
			
			if ($closest_ancestors{$site_index}){
				my $anc_node = $closest_ancestors{$site_index};
				#my $halfdepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()} - ($node->get_branch_length)/2; #19.09.2016 
				#my $fulldepth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()}; #19.09.2016 
			#	print " ancestor ".$anc_node->get_name(). " node ".$node->get_name()." depth $depth\n";
			#	push $static_subtree_info{$anc_node->get_name()}{$site_index}{"nodes"}, \$node;
				if (!$no_neighbour_changing || ($no_neighbour_changing && ! compare::is_neighbour_changing($self->{static_subs_on_node}{$node->get_name()}{$site_index}, 1))){
					if (!$no_leaves || ($no_leaves && !($node->is_terminal()))){
						push @{$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"synresearch"}}, $self->{static_subs_on_node}{$node->get_name()}{ $site_index}; 
						push @{$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"list"}}, $node->get_name();
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
   
   	sub lrt_visitor {
   		my $self = shift;
 		my $node = $_[0];

		if (!$node->is_root){
			my %closest_ancestors = %{$node->get_parent->get_generic("-closest_ancestors")};
			
			if (%closest_ancestors){
				foreach my $site_index(keys %closest_ancestors){ 
					my $anc_node = $closest_ancestors{$site_index};
					my $depth = $self->{static_distance_hash}{$anc_node->get_name()}{$node->get_name()};
					$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"lrt"}{$node->get_name()}[1] = $depth-($node->get_branch_length); #23.03 t_branch_start
					$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"lrt"}{$node->get_name()}[2] = $depth; #23.03 t_branch_end					
				}
				
				my @ancestors = keys %closest_ancestors;	
				foreach my $site_index(@ancestors){
					if (!($self->has_no_background_mutation($node, $site_index))){
						delete $closest_ancestors{$site_index};
					}
				}	
			}
			
			foreach my $site_index(keys %{$self->{static_subs_on_node}{$node->get_name()}}){
				if ($closest_ancestors{$site_index}){
					my $anc_node = $closest_ancestors{$site_index};
					$self->{static_subtree_info}{$anc_node->get_name()}{$site_index}{"lrt"}{$node->get_name()}[0] = 1;
				}
				$closest_ancestors{$site_index} = $node;
			}
			
			$node->set_generic("-closest_ancestors" => \%closest_ancestors);
		}

 	}
   
   sub predefined_groups_and_names {
 		my $self = shift;
 		my $prot = $self->{static_protein};
 		return Groups::get_predefined_groups_and_names_for_protein($prot, $self->mylength());
   }
   
   sub fake_predefined_groups_and_names {
   	 	my $self = shift;
   	 	my $exclude = shift;
 		my $prot = $self->{static_protein};
 		my $state = $self->{static_state};
 		return Groups::get_fake_predefined_groups_and_names_for_protein($prot, $self->mylength(), $exclude, $state);
   }

	sub protein_no_group {
 		my $self = shift;
 		my $prot = $self->{static_protein};
 		return Groups::get_no_groups_for_protein($prot, $self->mylength());
   }
   
 # changed at 21.09.2016  
   sub bin {
   	my $depth = $_[0];
   	my $step = $_[1];
   	
   	my $bin = int($depth/$step);
   	if ((int($depth/$step) == $depth/$step && $depth != 0) || $step == 1 || $step == 0.5){ # 0 goes to 0 bin, if step is 0.5 or 1, and to 1 bin otherwise
   		$bin -= 1;
   	}
   	return $bin+1;
   }
   

   

    

1;
