#!/usr/bin/perl

use Git::Repository;

package Codeversion;

my $version = undef;

sub get_version{
	if (! defined $version){
 		my $r = Git::Repository->new();
 		$version = $r->run( "rev-parse", "HEAD" );
	}
	return $version;
}
 1;