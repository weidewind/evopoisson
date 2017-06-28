#!/bin/bash

foldername="ksu_skip_stoppers"
perl poissonnolim.pl -p h3 --state nsyn --input little_ksu --step 0.00005 --skip_stoppers --restrictions 0.03,0.06 --simnumber 5000 --mutnum_control 0 --shuffler_type exp --maxmem 250000000 --output $foldername &>output/logs/ksudata_whole_skst
perl single_sites_nolim.pl -p h3 --state nsyn --input little_ksu --step 0.00005 --skip_stoppers --mutnum_control 0  --restrictions 0.03 --output $foldername &>output/logs/ksudata_sites_skst
		
perl grep_results.pl -i output/${foldername}/nsyn/maxpath_not_subtracted
perl grep_results.pl -i output/${foldername}/syn/maxpath_not_subtracted
exit 0
