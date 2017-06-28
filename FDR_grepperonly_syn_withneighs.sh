#!/bin/bash

foldername="little_ksu/fakes_littleksu_h3"


for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/${foldername}/${it}_fake/syn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -p h3 -r 0.03 -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl --state syn -i output/${foldername}
perl  metagrep_newsim_stats.pl --state syn -i output/${foldername}
exit 0
