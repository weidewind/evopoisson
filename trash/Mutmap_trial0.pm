#!/usr/bin/perl 

package Mutmap;

use strict;
use warnings;
use Cwd qw(abs_path cwd getcwd);
use File::Spec;
use Config::Simple;
use base qw(Class::Accessor); # inherits from Class::Accessor

use Bio::Phylo::IO;

# creates getters
Mutmap->mk_accessors(qw(_protein _tree _fasta _alignment_length _hash_of_nodes _state _subs_on_node _nodes_with_sub _background_subs_on_node _background_nodes_with_sub));


sub confManager {
	my $filepath;
	if ($_[0]){ $filepath = $_[0]; }
	else { $filepath = File::Spec->catfile(getcwd(), "data","default_config"); }
	my %config;
	if (-f $filepath)  { 
		Config::Simple->import_from($filepath, \%config);
		
		## checking config file consistency
		#die "only syn or nsyn can be used as the second argument; unknown $_[1] was used instead";
		 
		return %config;	 
	}
	else { die "Cannot open config file $filepath"; }
}


sub new {
    my $class = shift;
    my %config = confManager(shift);	
    my $tree = parse_tree(confManager->{tree});
    my @arr = parse_fasta(confManager->{fasta});
    my %fasta = %{$arr[0]};
    my $alignment_length = $arr[1];
    my $protein  = confManager->{protein};
    my %hash_of_nodes;    
        
    my @nodes = $tree -> get_nodes;
    foreach my $node(@nodes){
        #if ($node->is_root()) {next;}
        my $name = $node -> get_name();
        $hash_of_nodes{$name} = \$node;
    }
    
    my @mutmaps;
    my @bkg_mutmaps;
    my $state = confManager->{state};
    if($state eq "syn"){
        @mutmaps = synmutmap($tree, \%fasta);
        @bkg_mutmaps = codonmutmap($tree, \%fasta);

    } 
    elsif($state eq "nsyn"){
        @mutmaps = codonmutmap($tree, \%fasta);
        @bkg_mutmaps = synmutmap($tree, \%fasta);
    } 

    
    my $self = {
        _config => \%config,
        _tree => $tree,
        _fasta => \%fasta,
        _alignment_length => $alignment_length,
        _hash_of_nodes  => \%hash_of_nodes,
        _subs_on_node => \%{$mutmaps[0]},
        _nodes_with_sub => \%{$mutmaps[1]},
        _background_subs_on_node => \%{$bkg_mutmaps[0]},
        _background_nodes_with_sub => \%{$bkg_mutmaps[1]},
    };

    bless $self, $class;
    return $self;
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
		#my @nsyn = nsyn_substitutions($nodeseqs{$node->get_ancestors()->[0]->get_name()},
		#							  $nodeseqs{$name});
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


1;