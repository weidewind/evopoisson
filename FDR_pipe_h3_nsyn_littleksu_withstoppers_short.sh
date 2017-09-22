#!/bin/bash
called=`basename "$0"`
inputname="little_ksu"
prot="h3"
state="nsyn"
outputname="fakes_h3_withstoppers"
foldername="output/${inputname}/${outputname}"
mkdir -p $foldername
cp $called "${foldername}/${called}"

for i in {0..1}
do
		for number in {1..5}
		do
			it=$((number+5*i))
				sleep 3
				out="${outputname}/${it}_fake"
				( perl poissonnolim.pl -p $prot --state $state --input $inputname \
--step 0.00005 --restrictions 0,0.03,0.06 --fails_threshold 0.5 --mutnum_control 0 \
--skip_stoppers_in_simulation --simnumber 500  --fake --shuffler_type exp --maxmem 25000000 --onlysim \
--output $out &>output/logs/${outputname};
				perl groups_nolim.pl -p $prot --state $state --input $inputname \
--step 0.00005 --restrictions 0,0.03,0.06 --fails_threshold 0.5 --mutnum_control 0 \
--skip_stoppers_in_simulation \
--output $out &>"output/logs/${outputname}_groups" &
 				perl single_sites_nolim.pl -p $prot --state $state --input $inputname \
--step 0.00005 --restrictions 0 --fails_threshold 0.5 --mutnum_control 0 \
--skip_stoppers_in_simulation \
--output $out &>output/logs/${outputname}_sites) &
		done
		wait
		
done



for i in {0..1}
do
		for number in {1..5}
		do
			it=$((number+5*i))
				out="${foldername}/${it}_fake/${state}/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl -i ${foldername}
perl  metagrep_newsim_stats.pl -i ${foldername}
exit 0
