#!/bin/bash

while getopts p:d:n: option
do
		case "${option}" in
			p) PROTEIN=${OPTARG};;
			d) DIR=${OPTARG};;
			n) NUM=${OPTARG};;

		esac
done

called=`basename "$0"`
FDRFOLDER="output/${DIR}_${PROTEIN}_FDR"

for (( i=1; i<=$NUM; i++ ))
do
	OUT="${DIR}_${PROTEIN}_${i}"
	mv "output/${OUT}" "${FDRFOLDER}/${i}_fake"
done	

cd ..
for (( it=1; it<=$NUM; it++ ))
do
		out="${FDRFOLDER}/${it}_fake/nsyn/maxpath_not_subtracted/skip_stoppers"
		perl grep_newsim_stats.pl -i $out -r 0 &
		
done
wait

	
perl grep_fake_fdr.pl --skip_stoppers -r 0 -i ${FDRFOLDER}
perl  metagrep_newsim_stats.pl --skip_stoppers -r 0 -i ${FDRFOLDER}
exit 0