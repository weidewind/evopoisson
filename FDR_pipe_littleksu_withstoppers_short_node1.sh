#!/bin/bash
while getopts p:s:o: option
do
		case "${option}" in
			p) prot=${OPTARG};;
			s) state=${OPTARG};;
			o) omit=${OPTARG};;
		esac
done
called=`basename "$0"`
inputname="little_ksu"
outputname="fakes_${prot}_${state}_withstoppers"
foldername="output/${inputname}/${outputname}"
mkdir -p $foldername
cat $called >>"${foldername}/${called}"
echo "# protein ${prot} state ${state} omit ${omit}" >>"${foldername}/${called}"
for i in {0..25}
do
		for number in {1..4}
		do
			it=$((omit+number+4*i))
				sleep 3
				out="${outputname}/${it}_fake"
				( perl poissonnolim.pl -p $prot --state $state --input $inputname \
--step 0.00005 --restrictions 0,0.03,0.06 --mutnum_control 0 --fails_threshold 0.2 \
--simnumber 1000  --fake --shuffler_type exp --maxmem 5000000 --onlysim \
--output $out &>output/logs/${outputname};
				perl groups_nolim.pl -p $prot --state $state --input $inputname \
--step 0.00005 --restrictions 0,0.03,0.06 --mutnum_control 0 --fails_threshold 0.2 \
--output $out &>"output/logs/${outputname}_groups" &
 				perl single_sites_nolim.pl -p $prot --state $state --input $inputname \
--step 0.00005 --restrictions 0 --mutnum_control 0 --fails_threshold 0.2 \
--output $out &>output/logs/${outputname}_sites) &
		done
		wait
		
done



for i in {0..14}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="${foldername}/${it}_fake/${state}/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out -p $prot --state $state -r 0 &
				
		done
		wait
done
	
perl grep_fake_fdr.pl -i ${foldername} --state $state
perl  metagrep_newsim_stats.pl -i ${foldername} --state $state
exit 0