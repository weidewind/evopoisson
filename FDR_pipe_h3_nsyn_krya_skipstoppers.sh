#!/bin/bash

foldername="fakes_h3_krya_skipstoppers"
for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				sleep 3
				out="${foldername}/${it}_fake"
				( perl poissonnolim.pl -p h3 --state nsyn --step 0.5 --restrictions 50,100 --onlysim --simnumber 1000 --skip_stoppers --mutnum_control 0 --fake --shuffler_type exp --maxmem 25000000 --output $out &>output/logs/fakes_h3_nsyn_wst;
				  perl groups_nolim.pl -p h3 --state nsyn --restrictions 50,100 --mutnum_control 0 --skip_stoppers --output $out &>output/logs/fakes_h3_nsyn_groups_wst &
 				perl single_sites_nolim.pl -p h3 --state nsyn --step 0.5 --restrictions 50 --skip_stoppers --mutnum_control 0 --output $out &>output/logs/fakes_h3_nsyn_sites_wst) &
		done
		wait
		
done



for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/${foldername}/${it}_fake/nsyn/maxpath_not_subtracted/skip_stoppers"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl --skip_stoppers -i output/little_ksu/${foldername}
perl  metagrep_newsim_stats.pl --skip_stoppers -i output/little_ksu/${foldername}
exit 0
