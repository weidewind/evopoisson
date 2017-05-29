#!/bin/bash

foldername="nolim_massive_fakes_n1_after_serverbugfix"
for i in {1..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				sleep 3
				out="${foldername}/${it}_fake"
				perl poissonnolim.pl -p n1 --state syn --no_neighbour_changing --simnumber 1000 --mutnum_control 0 --fake --shuffler_type exp --maxmem 50000000 --output $out &
		done
		wait
		
		for number in {1..10}
		do
			it=$((number+10*i))
				out="${foldername}/${it}_fake"
				perl single_sites_nolim.pl -p n1 --state syn --no_neighbour_changing --mutnum_control 0 --output $out&
		done
		wait
done



for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/${foldername}/${it}_fake/syn/maxpath_not_subtracted/no_neighbour_changing"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl --state syn --no_neighbour_changing -i output/${foldername}/syn
perl  metagrep_newsim_stats.pl --state syn --no_neighbour_changing -i output/${foldername}/syn
exit 0
