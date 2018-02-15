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
outputname="realdata_${prot}_${state}_withstoppers"
foldername="output/${inputname}/${outputname}"
mkdir -p $foldername
cat $called >>"${foldername}/${called}"
echo "# protein ${prot} state ${state}" >>"${foldername}/${called}"

perl poissonnolim.pl -p $prot --state $state --input $inputname --step 0.00005 --restrictions 0,0.03,0.06 \
--onlysim --fails_threshold 0.2 --mutnum_control 0 \
--simnumber 10000 --shuffler_type exp --maxmem 25000000 --output $outputname &>output/logs/${outputname};
perl groups_nolim.pl -p $prot --state $state --input $inputname  --step 0.00005 --restrictions 0,0.03,0.06 \
--mutnum_control 0 --fails_threshold 0.2 --output $outputname &>"output/logs/${outputname}_groups" &
perl single_sites_nolim.pl -p $prot --state $state --input $inputname --step 0.00005 --restrictions 0 \
--mutnum_control 0 --fails_threshold 0.2 --output $outputname &>output/logs/${outputname}_sites &
wait

