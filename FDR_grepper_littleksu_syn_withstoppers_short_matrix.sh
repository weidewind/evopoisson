#!/bin/bash
while getopts p: option
do
		case "${option}" in
			p) prot=${OPTARG};;
		esac
done
called=`basename "$0"`
inputname="little_ksu"
state="syn"
outputname="fakes_matrix_${prot}_syn_withstoppers"
foldername="output/${inputname}/${outputname}"
mkdir -p $foldername
cat $called >>"${foldername}/${called}"
echo "# protein ${prot} state syn" >>"${foldername}/${called}"



for i in {0..20}
do
		for number in {1..5}
		do
			it=$((number+5*i))
				out="${foldername}/${it}_fake/${state}/maxpath_not_subtracted"
				perl grep_newsim_stats.pl -i $out -p $prot --state $state -r 0 &>output/logs/${outputname}_grepper &
				
		done
		wait
done
	
perl grep_fake_fdr.pl -i ${foldername} --state $state
perl  metagrep_newsim_stats.pl -i ${foldername} --state $state
exit 0
