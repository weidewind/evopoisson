#!/bin/bash
for i in {0..49}
do
	for number in {1..10}
	do
			it=$((number+20*i))
			perl FDR_all.pl -p h1 --state nsyn --output fake_no_restriction --num $it --restriction 0 >output/fake_new_clones/nsyn/maxpath_not_subtracted/logs &
	done
	wait
done	
exit 0
