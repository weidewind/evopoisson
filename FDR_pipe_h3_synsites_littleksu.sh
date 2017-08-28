#!/bin/bash

foldername="fakes_littleksu_h3"
for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				sleep 3
				out="${foldername}/${it}_fake"
 				perl single_sites_nolim.pl -p h3 --state syn --input little_ksu --step 0.00005 --restrictions 0.03 --mutnum_control 0 --output $out &>output/logs/fakes_littleksu_synsites &
		done
		wait
		
done



for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/little_ksu/${foldername}/${it}_fake/syn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl -i output/little_ksu/${foldername}
perl  metagrep_newsim_stats.pl -i output/little_ksu/${foldername}
exit 0
