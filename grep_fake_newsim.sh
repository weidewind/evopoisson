#!/bin/bash

for i in {0..5}
do
		for number in {1..10}
		do
			it=$((number+10*i))
			if (($it < 111)); then
				out="output/onestrip_debugged/${it}_fake/nsyn/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out &
			fi	
		done
		wait
done
	
exit 0