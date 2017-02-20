#!/bin/bash

foldername="exp_shuffler"
for i in {0..5}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="${foldername}/${it}_fake"
				perl poisson.pl -p h1 --state nsyn --simnumber 100 --mutnum_control 0 --fake --maxmem 8000000 --output $out >>output/fakeslog_${foldername} &
		done
		wait
		
		for number in {1..10}
		do
			it=$((number+10*i))
				out="${foldername}/${it}_fake"
				perl single_sites.pl -p h1 --state nsyn --mutnum_control 0 --output $out >>output/fakesiteslog_${foldername} &
		done
		wait
done



for i in {0..5}
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
