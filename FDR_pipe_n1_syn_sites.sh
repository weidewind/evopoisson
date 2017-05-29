#!/bin/bash

foldername="nolim_massive_fakes_n1"
for i in {0..10}
do
		for number in {1..15}
		do
			it=$((number+15*i+210))
				out="${foldername}/${it}_fake"
				perl single_sites_nolim.pl -p n1 --state syn --no_neighbour_changing --mutnum_control 0 --output $out&
		done
		wait
done



for i in {0..10}
do
		for number in {1..15}
		do
			it=$((number+15*i+210))
				out="output/${foldername}/${it}_fake/syn/maxpath_not_subtracted/no_neighbour_changing"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl --state syn --no_neighbour_changing -i output/${foldername}/syn
perl  metagrep_newsim_stats.pl --state syn --no_neighbour_changing -i output/${foldername}/syn
exit 0
