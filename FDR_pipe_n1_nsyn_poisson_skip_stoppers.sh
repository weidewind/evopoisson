#!/bin/bash

foldername="nolim_massive_fakes_skipstoppers_poisson_n1"
for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="${foldername}/${it}_fake"
				sleep 3
				( perl poissonnolim.pl -p n1 --state nsyn --skip_stoppers --distrpoisson --simnumber 1000 --mutnum_control 0 --fake --shuffler_type exp --maxmem 60000000 --output $out ;
				perl single_sites_nolim.pl -p n1 --state nsyn --mutnum_control 0 --output $out ) &
		done
		wait
done



for i in {0..10}
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
