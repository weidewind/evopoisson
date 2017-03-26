#!/bin/bash

foldername="fake_exp_shuffler_after_middle_debugging"
for i in {0..3}
do	
		for number in {1..10}
		do
			it=$((number+10*i))
				out="${foldername}/${it}_fake"
				perl single_sites.pl -p h1 --state nsyn --mutnum_control 0 --output $out >>output/fakesiteslog_${foldername} &
		done
		wait
done



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

perl  metagrep_newsim_stats.pl -i output/${foldername}
exit 0
