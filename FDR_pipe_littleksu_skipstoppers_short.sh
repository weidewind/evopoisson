#!/bin/bash
while getopts p:s: option
do
		case "${option}" in
			p) prot=${OPTARG};;
			s) state=${OPTARG};;
		esac
done
called=`basename "$0"`
inputname="little_ksu"
outputname="fakes_${prot}_skipstoppers"
foldername="output/${inputname}/${outputname}"
mkdir -p $foldername
cat $called >>"${foldername}/${called}"
echo "# protein ${prot} state ${state}" >>"${foldername}/${called}"
for i in {0..14}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				sleep 3
				out="${outputname}/${it}_fake"
				( perl poissonnolim.pl -p $prot --state $state --input $inputname \
--skip_stoppers --step 0.00005 --restrictions 0,0.03,0.06 --mutnum_control 0 --fails_threshold 0.2 \
--skip_stoppers_in_simulation --simnumber 1000  --fake --shuffler_type exp --maxmem 25000000 --onlysim \
--output $out &>output/logs/${outputname};
				perl groups_nolim.pl -p $prot --state $state --input $inputname \
--skip_stoppers --step 0.00005 --restrictions 0,0.03,0.06 --mutnum_control 0 --fails_threshold 0.2 \
--skip_stoppers_in_simulation \
--output $out &>"output/logs/${outputname}_groups" &
 				perl single_sites_nolim.pl -p $prot --state $state --input $inputname \
--skip_stoppers --step 0.00005 --restrictions 0 --mutnum_control 0 --fails_threshold 0.2 \
--skip_stoppers_in_simulation \
--output $out &>output/logs/${outputname}_sites) &
		done
		wait
		
done



for i in {0..14}
do
		for number in {1..10}
		do
			it=$((number+10*i))
				out="${foldername}/${it}_fake/${state}/maxpath_not_subtracted/skip_stoppers"
				perl grep_newsim_stats.pl -i $out -p $prot --state $state -r 0 &
				
		done
		wait
done
	
perl grep_fake_fdr.pl -i ${foldername} --skip_stoppers --state $state
perl  metagrep_newsim_stats.pl -i ${foldername} --skip_stoppers --state $state
exit 0
