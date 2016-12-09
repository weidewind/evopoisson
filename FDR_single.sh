#!/bin/bash
for i in {0..29}
do
	for number in {1..15}
	do
			it=$((number+20*i))
			perl FDR_single_sites.pl -p h1 --state nsyn --output fake_new_iterations --num $it  >output/fake_new_iterations/nsyn/maxpath_not_subtracted/logs_single &
	done
	wait
done	
exit 0
