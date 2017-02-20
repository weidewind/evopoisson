#!/bin/bash

for i in {0..5}
do
		for number in {1..10}
		do
			it=$((number+10*i))
			if (($it < 111)); then
				out="onestrip_no_restriction/${it}_fake"
				perl single_sites.pl -p h1 --state nsyn --mutnum_control 0 --output $out >>output/fakesiteslog_onestrip_norestriction &
			fi	
		done
		wait
done

for i in {0..5}
do
		for number in {1..10}
		do
			it=$((number+10*i))
			if (($it < 111)); then
				out="output/onestrip_no_restriction/${it}_fake/nsyn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out &
			fi	
		done
		wait
done
	
exit 0
