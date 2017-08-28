#!/bin/bash
called =`basename "$0"`
foldername="real_h1_krya_skipstoppers_all"
cp $called "output/${foldername}/bashscript"
perl poissonnolim.pl -p h1 --state nsyn --step 0.5 --restrictions 0 --onlysim --simnumber 10000 --skip_stoppers --mutnum_control 0 --shuffler_type exp --maxmem 25000000 --output $foldername &>output/logs/h1_nsyn_all
perl groups_nolim.pl -p h1 --state nsyn --step 0.5 --restrictions 0 --mutnum_control 0 --skip_stoppers --output $foldername &>output/logs/h1_nsyn_groups_all &
perl single_sites_nolim.pl -p h1 --state nsyn --step 0.5 --restrictions 0 --skip_stoppers --mutnum_control 0 --output $foldername &>output/logs/h1_nsyn_sites_all
