#!/usr/bin/perl
## This class provides a map of mutations (node_name -> array of sites that mutated on the branch),  
## given a tree and sequences of all nodes (including internal ones)

## 

package Locker;


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

sub try_and_start {
	my $self = shift;
	my $start = 0;
	my $lockfile = $self->{lockfile};
	my $maxmem = $self->{maxmem};
	open LOCK, "+<$lockfile" or die; # open for update (you can read and update!)
      	 flock (LOCK, 2); # it also checks for any locks and waits if finds any
      	 print "start locked\n";
      	 my $proc_memusage = <LOCK>; # undef on start and while first process is running 
      	 chomp ($proc_memusage);
      	#unless ($proc_memusage && $proc_memusage ne '' && $proc_memusage ne "\n" )  {$proc_memusage = 0;}
      	 my $num_running = <LOCK>; # undef on start
      	 chomp($num_running);
      	#unless ($num_running && $num_running ne '' && $num_running ne "\n") {$num_running = 0;}
      	 if (!$num_running || ($proc_memusage && $maxmem > ($num_running+1)*$proc_memusage)){
      	 	seek LOCK, 0, SEEK_SET;
      	 	truncate(LOCK,0); # empty the file
      	 	print LOCK $proc_memusage."\n".($num_running+1);
      	 	if (!$num_running){print "!$num_running is true\n";}
      	 	if (($proc_memusage && $maxmem > ($num_running+1)*$proc_memusage)){print "eq is true\n";}
      	 	print "maxmem  $maxmem > ".($num_running+1)*$proc_memusage." where memusage = $proc_memusage \n";
      	 	print "Starting process: already printed to lock ".$proc_memusage." ".($num_running+1)."eof\n";
      	 	close LOCK;
      	 	print "start unlocked\n";
         	$start = 1;
         	
      	 }
      	 else {
      	 	close LOCK;
      	 	print "start unlocked\n";
      	 	print "process not ready, going to sleep..\n";
      	 	sleep (15);
      	 }
      	 return $start;
}

sub delete_entry {
		my $self = shift;
		my $lockfile = $self->{lockfile};
	    open LOCK, "+<$lockfile" or die; # open for update
      	flock (LOCK, 2); # it also checks for any locks and waits if finds any
      	print "fin locked\n";
      	my $proc_memusage = <LOCK>; # undef on start and while first process is running 
	    chomp ($proc_memusage);
	    my $num_running = <LOCK>; # undef on start
	    chomp($num_running);
      	#unless ($num_running && $num_running ne '' && $num_running ne "\n") {$num_running = 0;}
      	seek LOCK, 0, SEEK_SET;
      	truncate(LOCK,0);
      	print LOCK $proc_memusage."\n".($num_running-1);
      	print "Exiting: already printed to lock ".$proc_memusage." ".($num_running-1)."eof\n";
      	close LOCK;
      	print "fin unlocked\n";
}

sub print_memusage {
		my $self = shift;
		my $lockfile = $self->{lockfile};
	    open LOCK, "+<$lockfile" or die; # open for update
      	flock (LOCK, 2); # it also checks for any locks and waits if finds any
      	print "pm locked\n";
      	my $proc_memusage = <LOCK>; # undef on start and while first process is running 
	    chomp ($proc_memusage);
	    my $num_running = <LOCK>; # undef on start
	    chomp($num_running);
	    seek LOCK, 0, SEEK_SET;
      	truncate(LOCK,0);
      	my $mu = Memory::Usage->new();
		$mu->record('inside iterations_gulp');
		my $proc_memusage = $mu->state()->[0]->[2]; # virtual memory size
      	print LOCK $proc_memusage."\n".$num_running;
      	print "Iterating: already printed to lock ".$proc_memusage." ".$num_running."eof\n";
      	close LOCK;
      	print "pm unlocked\n";
}

1;