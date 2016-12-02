#!/bin/bash
for i in {0..25}
do
	for number in {1..20}
	do
			it=$((number+20*i))
			perl FDR_all.pl -p h1 --state nsyn --output fake_new_clones --num $it >output/fake_new_clones/nsyn/maxpath_not_subtracted/logs &
	done
	wait
done	
exit 0
