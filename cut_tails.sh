#!/bin/bash
for i in {1..33}
do 
	head -n $(grep -n '^>' h1_for_enrichment_$i | tail -n 1 | cut -d':' -f1) h1_for_enrichment_$i >temp
	head -n -1 temp >h1_for_enrichment_$i 
done
