#!/bin/bash

foldername="little_test_binback"

perl single_sites_nolim.pl -p h3 --state nsyn --step 0.0001 --input little_ksu --mutnum_control 0  --restrictions 0.03 --output $foldername &>output/logs/ksudata_sites &
perl single_sites_nolim.pl -p h3 --state syn --step 0.0001 --input little_ksu --restrictions 0.03 --mutnum_control 0 --output $foldername &>output/logs/ksudata_sites_syn
wait		
perl grep_results.pl -i output/${foldername}/nsyn/maxpath_not_subtracted
perl grep_results.pl -i output/${foldername}/syn/maxpath_not_subtracted
exit 0
