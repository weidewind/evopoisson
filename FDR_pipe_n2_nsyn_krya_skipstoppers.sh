#!/bin/bash

foldername="fakes_n2_nsyn_krya_skipstoppers"
for i in {0..9}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				sleep 3
				out="${foldername}/${it}_fake"
				( perl poissonnolim.pl -p n2 --state nsyn --step 0.5 --restrictions 50,100 --onlysim --simnumber 1000 --mutnum_control 0 --fake --shuffler_type exp --maxmem 3500000 --output $out &>output/logs/fakes_n2_nsyn_wst;
				  perl groups_nolim.pl -p n2 --state nsyn --restrictions 50,100 --mutnum_control 0 --output $out &>output/logs/fakes_n2_nsyn_groups_wst &
 				perl single_sites_nolim.pl -p n2 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output $out &>output/logs/fakes_n2_nsyn_sites_wst) &
		done
		wait
		
done



for i in {0..9}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="output/${foldername}/${it}_fake/nsyn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl -i output/${foldername}
perl  metagrep_newsim_stats.pl -i output/${foldername}
exit 0