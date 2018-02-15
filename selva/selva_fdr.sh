#!/bin/bash

while getopts c:p:d:n:u option
do
		case "${option}" in
			c) CONFIG=${OPTARG};;
			p) PROTEIN=${OPTARG};;
			d) DIR=${OPTARG};;
			n) NUM=${OPTARG};;
			u) UPDATE_TREE=true;;
		esac
done

called=`basename "$0"`
FDRFOLDER="output/${DIR}_${PROTEIN}_FDR"
cd ..
mkdir "${FDRFOLDER}"
cp "selva/${called}" "${FDRFOLDER}/${called}"
echo "# selva_fdr options: CONFIG=${CONFIG} PROTEIN=${PROTEIN} DIR=${FDRFOLDER} NUM=${NUM} UPDATE_TREE=${UPDATE_TREE}\n" >>"${FDRFOLDER}/${called}"
cd selva

for (( i=1; i<=$NUM; i++ ))
do
	OUT="${DIR}_${PROTEIN}_${i}"
	if [ $UPDATE_TREE ] ; then
		./selva_pipe.sh -p $PROTEIN -c $CONFIG -d $OUT -u
	else
		./selva_pipe.sh -p $PROTEIN -c $CONFIG -d $OUT
	fi	
	( cd .. ;
	perl poissonnolim.pl -p $PROTEIN --state nsyn --input $OUT \
--step 0.5 --restrictions 0,3,6 --mutnum_control 0 \
--skip_stoppers_in_simulation --skip_stoppers --simnumber 1000 --shuffler_type exp --maxmem 1500000 --onlysim \
&>output/logs/${OUT};
	perl groups_nolim.pl -p $PROTEIN --state nsyn --input $OUT \
--step 0.5 --restrictions 0,3,6 --mutnum_control 0 --stat_types mean,median \
--skip_stoppers_in_simulation --skip_stoppers \
&>"output/logs/${OUT}_groups" &
 	perl single_sites_nolim.pl -p $PROTEIN --state nsyn --input $OUT \
--step 0.5 --restrictions 0 --mutnum_control 0 \
--skip_stoppers_in_simulation --skip_stoppers \
&>output/logs/${OUT}_sites &
	wait ) &
done
wait

#move each fake into fdrfolder and change their names
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
	
	

