#!/bin/bash

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