#!/bin/bash

foldername="little_ksudata"
perl poissonnolim.pl -p h3 --state syn --input little_ksu  --no_neighbour_changing --step 0.00005 --restrictions 0.03,0.06,0.09 --simnumber 10000 --mutnum_control 0 --shuffler_type exp --maxmem 250000000 --output $foldername &>output/logs/ksudata_whole_syn_noneigh
perl single_sites_nolim.pl -p h3 --state syn --input little_ksu --no_neighbour_changing --step 0.00005 --restrictions 0.03 --mutnum_control 0 --output $foldername &>output/logs/ksudata_sites_noneigh_syn

perl grep_results.pl -i output/${foldername}/syn/maxpath_not_subtracted/no_neighbour_changing
exit 0
