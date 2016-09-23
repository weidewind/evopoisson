#!/usr/bin/perl

use Git::Repository;

package Codeversion;

sub get_version{
 	my $r = Git::Repository->new();
 	my $output = $r->run( "rev-parse", "HEAD" );
 	return $output;
}
 1;