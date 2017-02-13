#!/bin/bash

for i in {0..10}
do
		for number in {1..10}
		do
			it=$((number+10*i))
			out="${it}_fake"
			perl poisson.pl -p h1 --state nsyn --simnumber 300 --mutnum_control 0 --fake --maxmem 6000000 --output $out >>output/fakeslog &
		done
		wait
done	
exit 0