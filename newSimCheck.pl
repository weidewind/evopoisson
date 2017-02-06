#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmap;


my $args = { bigtag => "shuffle_fake", protein => "h1", state => "nsyn", subtract_tallest => 0,    mutnum_control => 0, fromfile => 1}; 
$mutmap = Mutmap->new($args);
if (exists $mutmap->{realdata}{subtree_info}){
	print "yeah\n";
}
open GULP, "</export/home/popova/workspace/evopoisson/output/shuffle_fake/nsyn/maxpath_not_subtracted/h1/h1_for_enrichment_44";

while (<GULP>){
	if ($_ =~ /^s/ ){
		my ($a, $site, $b, $nodename, $c, $maxdepth, $d, $mutnum)= split(/\s+/, $_);
		print $site."__".$nodename."\n";
		my $rm = $mutmap->{realdata}{subtree_info}{$nodename}{$site}{"maxdepth"};
		if (exists $mutmap->{realdata}{subtree_info}{$nodename}{$site}){
			foreach my $key(keys %{$mutmap->{realdata}{static_subtree_info}{$nodename}{$site}}){
				print $key."\t";
			}
			print "\n";
		}
		print "real maxdepth ".$rm;
		print " sim maxdepth ".$maxdepth."\n";
		if ($maxdepth > $rm){
			print "Ooops!\n";
		}
		my $rn = $mutmap->{realdata}{subtree_info}{$nodename}{$site}{"totmuts"};
		print "real mutnum ".$rn;
		print " sim mutnum ".$mutnum."\n";
		if ($mutnum != $rn){
			print "Ouch!\n";
		}
	}

}
close GULP;
