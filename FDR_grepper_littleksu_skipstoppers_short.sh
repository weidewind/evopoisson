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

for i in {0..9}
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
