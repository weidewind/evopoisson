#!/bin/bash
for i in {3..6}
do
	for number in {1..20}
	do
			it=$((number+20*i))
			perl FDR_single_sites.pl -p h1 --state nsyn --output fake_new_iterations --num $it --mutnum_control >output/fake_new_iterations/nsyn/maxpath_not_subtracted/logs_single_mutnum &
	done
	wait
done	
exit 0
