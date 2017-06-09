#!/bin/bash

foldername="ksudata"
(perl poissonnolim.pl -p h3 --state nsyn --input ksu --restrictions 0.03,0.06,0.09 --simnumber 2000 --mutnum_control 0 --shuffler_type exp --maxmem 250000000 --output $foldername &>output/logs/ksudata_whole ;
perl single_sites_nolim.pl -p h3 --state nsyn input ksu --mutnum_control 0  --restrictions 0.03,0.06,0.09 --output $foldername &>output/logs/ksudata_sites ) &
perl poissonnolim.pl -p h3 --state syn --input ksu  --restrictions 0.03,0.06,0.09 --simnumber 2000 --mutnum_control 0 --shuffler_type exp --maxmem 250000000 --output $foldername &>output/logs/ksudata_whole_syn
perl single_sites_nolim.pl -p h3 --state syn input ksu --restrictions 0.03,0.06,0.09 --mutnum_control 0 --output $foldername &>output/logs/ksudata_sites_syn
		
perl grep_results.pl -i output/${foldername}/nsyn/maxpath_not_subtracted
perl grep_results.pl -i output/${foldername}/syn/maxpath_not_subtracted
exit 0
