#!/bin/bash

foldername="fake_onestrip_after_middle"
for i in {0..3}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/${foldername}/${it}_fake/nsyn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_FDR.pl -i output/${foldername}
perl  metagrep_newsim_stats.pl -i output/${foldername}
exit 0