#!/bin/bash

foldername="nolim_massive_fakes_n2"


for i in {0..20}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/${foldername}/${it}_fake/syn/maxpath_not_subtracted/no_neighbour_changing"
				perl grep_newsim_stats.pl -p n2 -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl --state syn --no_neighbour_changing -i output/${foldername}
perl  metagrep_newsim_stats.pl --state syn --no_neighbour_changing -i output/${foldername}
exit 0
