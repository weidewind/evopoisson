#!/bin/bash
called=`basename "$0"`
inputname="little_ksu"
prot="h3"
state="nsyn"
outputname="real_skipstoppers"
foldername="output/${outputname}"
mkdir -p $foldername
cp $called "${foldername}/${called}"

perl poissonnolim.pl -p $prot --state $state --input $inputname --step 0.00005 --restrictions 0,0.03,0.06 \
--skip_stoppers --skip_stoppers_in_simulation --onlysim --mutnum_control 0 \
--simnumber 10000 --shuffler_type exp --maxmem 25000000 --output $outputname &>output/logs/${outputname};
perl groups_nolim.pl -p $prot --state $state --input $inputname  --step 0.00005 --restrictions 0,0.03,0.06 \
--skip_stoppers --mutnum_control 0 --fails_threshold 0.5 --output $outputname &>"output/logs/${outputname}_groups" &
perl single_sites_nolim.pl -p $prot --state $state --input $inputname --step 0.00005 --restrictions 0 \
--skip_stoppers --mutnum_control 0 --output $outputname &>output/logs/${outputname}_sites &
wait

