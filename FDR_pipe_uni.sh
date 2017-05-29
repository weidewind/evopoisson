#!/bin/bash

foldername="nolim_massive_fakes_h3_after_serverbugfix"
protein="h3"
state="syn"
total_fakes=100
procnum=10
start_outer_index=0
end_outer_index=start_index+total_fakes/procnum
for i in {$start_outer_index..$end_outer_index}
do
		for number in {1..$procnum}
		do
			it=$((number+procnum*i))
				sleep 3
				out="${foldername}/${it}_fake"
				( perl poissonnolim.pl -p $protein --state $state --no_neighbour_changing --simnumber 1000 --mutnum_control 0 --fake --shuffler_type exp --maxmem 10000000 --output $out ;
				perl single_sites_nolim.pl -p $protein --state $state --no_neighbour_changing --mutnum_control 0 --output $out ) &
		done
		wait
done



for i in {$start_outer_index..$end_outer_index}
do
		for number in {1..procnum}
		do
			it=$((number+procnum*i))
				out="output/${foldername}/${it}_fake/syn/maxpath_not_subtracted/no_neighbour_changing"
				perl grep_newsim_stats.pl -i $out &
				
		done
		wait
done
	
perl grep_fake_fdr.pl --state $state --no_neighbour_changing -i output/${foldername}/${state}
perl  metagrep_newsim_stats.pl --state $state --no_neighbour_changing -i output/${foldername}/${state}
exit 0
