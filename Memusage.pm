#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

## 

package Memusage;


use File::Path qw(make_path remove_tree);
use Memory::Usage;



sub new {
		my ($class, $mutmap, $maxmem) = @_;
		my $lockfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}, "lock");
		my $self = {lockfile =>  $lockfile, maxmem => $maxmem};
		open LOCK, ">$lockfile" or die $!;
		close LOCK;	
		bless $self, $class;
		return $self;
}

# used inside iterations_gulp, does not create new lockfile
sub get_locker {
	my ($class, $mutmap) = @_;
		my $lockfile = File::Spec->catfile($mutmap->{static_output_base}, $mutmap->{static_protein}, "lock");
		my $self = {lockfile =>  $lockfile};
		bless $self, $class;
		return $self;
}

sub print_memusage {
		my $self = shift;
		my $lockfile = $self->{lockfile};
	    open MEM, ">$lockfile" or die $!; 
      	my $mu = Memory::Usage->new();
		$mu->record('inside iterations_gulp');
		my $proc_memusage = $mu->state()->[0]->[2]; # virtual memory size
      	print MEM $proc_memusage."\n";
      	print "Proc memusage is ".$proc_memusage."\n";
      	close MEM;
}

sub get_memusage {
		my $self = shift;
		my $lockfile = $self->{lockfile};
	    open MEM, "<$lockfile" or die $!." $lockfile"; 
      	my $memusage = <MEM>;
      	close MEM;
      	return $memusage;
}

1;