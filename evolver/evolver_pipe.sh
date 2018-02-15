#!/bin/bash

while getopts d:n: option
do
		case "${option}" in
#			p) PROTEIN=${OPTARG};;
			d) DIR=${OPTARG};;
			n) NUM=${OPTARG};;
		esac
done

INPUT="little_ksu"
PROTEIN="h1"

called=`basename "$0"`
FDRFOLDER="output/${DIR}_FDR"
mkdir "${FDRFOLDER}"
cp "evolver/${called}" "${FDRFOLDER}/${called}"
echo "# evolver_fdr options: PROTEIN=${PROTEIN} DIR=${FDRFOLDER} NUM=${NUM}\n" >>"${FDRFOLDER}/${called}"

perl EvolverParser.pl -d "evolver/${DIR}" -t "../data/${INPUT}/${PROTEIN}.l.r.newick"

for (( i=1; i<=$NUM; i++ ))
do 
	DATDIR="data/${DIR}_${i}"
	mkdir $DATDIR
	cp "evolver/${DIR}/evolver_${i}.fasta" "${DATDIR}/${PROTEIN}.all.fa"
	cp "data/${INPUT}/${PROTEIN}.l.r.newick" "${DATDIR}/${PROTEIN}.l.r.newick"
	if [ $UPDATE_TREE ]
	then 
		echo "Branch lengths will be updated\n"
		perl print_nsyn_tree.pl -p $PROTEIN --input "${DIR}_${i}"
		mv "${DATDIR}/${PROTEIN}.l.r.updated.newick" "${DATDIR}/${PROTEIN}.l.r.newick"
	fi
done

for (( i=1; i<=$NUM; i++ ))
do
	OUT="${DIR}_${i}"
(
	perl poissonnolim.pl -p $PROTEIN --state nsyn --input $OUT \
--step 0.00005 --restrictions 0,0.03,0.06 --mutnum_control 0 \
--skip_stoppers_in_simulation --skip_stoppers --simnumber 1000 --shuffler_type exp --maxmem 1500000 --onlysim \
&>output/logs/${OUT};
	perl groups_nolim.pl -p $PROTEIN --state nsyn --input $OUT \
--step 0.00005 --restrictions 0,0.03,0.06 --mutnum_control 0 --stat_types mean,median \
--skip_stoppers_in_simulation --skip_stoppers \
&>"output/logs/${OUT}_groups" &
 	perl single_sites_nolim.pl -p $PROTEIN --state nsyn --input $OUT \
--step 0.00005 --restrictions 0 --mutnum_control 0 \
--skip_stoppers_in_simulation --skip_stoppers \
&>output/logs/${OUT}_sites &
	wait ) &
done
wait

#move each fake into fdrfolder and change their names
for (( i=1; i<=$NUM; i++ ))
do
	OUT="${DIR}_${i}"
	mv "output/${OUT}" "${FDRFOLDER}/${i}_fake"
done	

for (( it=1; it<=$NUM; it++ ))
do
		out="${FDRFOLDER}/${it}_fake/nsyn/maxpath_not_subtracted/skip_stoppers"
		perl grep_newsim_stats.pl -i $out -r 0 &
		
done
wait

	
perl grep_fake_fdr.pl --skip_stoppers -r 0 -i ${FDRFOLDER}
perl  metagrep_newsim_stats.pl --skip_stoppers -r 0 -i ${FDRFOLDER}
exit 0
