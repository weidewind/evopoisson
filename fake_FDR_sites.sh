#!/bin/bash

for i in {0..3}
do
		for number in {1..35}
		do
			it=$((number+35*i))
			if (($it < 111)); then
				out="${it}_fake"
				perl single_sites.pl -p h1 --state nsyn --mutnum_control 0 --output $out >>output/fakesiteslog &
			fi	
		done
		wait
done

for i in {0..3}
do
		for number in {1..35}
		do
			it=$((number+35*i))
			if (($it < 111)); then
				out="output/${it}_fake/nsyn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out &
			fi	
		done
		wait
done
	
exit 0
