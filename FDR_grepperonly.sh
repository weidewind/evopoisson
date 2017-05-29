#!/bin/bash

foldername="nolim_massive_fakes_n1"


for i in {0..20}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/${foldername}/${it}_fake/nsyn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -p n1 -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl -i output/${foldername}
perl  metagrep_newsim_stats.pl -i output/${foldername}
exit 0
