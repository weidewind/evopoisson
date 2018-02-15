#!/bin/bash
FDRFOLDER="output/nabifdr_h1_FDR"
NUM=60
cd ..
for (( it=1; it<=$NUM; it++ ))
do
		out="${FDRFOLDER}/${it}_fake/nsyn/maxpath_not_subtracted/skip_stoppers"
		perl grep_newsim_stats.pl -i $out -r 0 &
		
done
wait

	
perl grep_fake_fdr.pl --skip_stoppers -i ${FDRFOLDER} -r 0
perl  metagrep_newsim_stats.pl --skip_stoppers -i ${FDRFOLDER} -r 0 
exit 0