#!/bin/bash
for i in {0..49}
do
	for number in {1..10}
	do
			it=$((number+20*i))
			perl FDR_all.pl -p h1 --state nsyn --output shuffle_test --num $it --restriction 50 >output/shuffle_test/nsyn/maxpath_not_subtracted/logs &
	done
	wait
done	
exit 0
