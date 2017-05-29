#!/bin/bash

foldername="nolim_massive_fakes_h3"
for i in {11..21}
do
		for number in {1..5}
		do
			it=$((number+5*i))
				out="${foldername}/${it}_fake"
				perl poissonnolim.pl -p h3 --state syn --no_neighbour_changing --simnumber 1000 --mutnum_control 0 --fake --shuffler_type exp --maxmem 40000000 --output $out &
		done
		wait
		
		for number in {1..5}
		do
			it=$((number+5*i))
				out="${foldername}/${it}_fake"
				perl single_sites_nolim.pl -p h3 --state syn --no_neighbour_changing --mutnum_control 0 --output $out&
		done
		wait
done



for i in {11..21}
do
		for number in {1..5}
		do
			it=$((number+4*i))
				out="output/${foldername}/${it}_fake/syn/maxpath_not_subtracted/no_neighbour_changing"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl --state syn --no_neighbour_changing -i output/${foldername}/syn
perl  metagrep_newsim_stats.pl --state syn --no_neighbour_changing -i output/${foldername}/syn
exit 0
