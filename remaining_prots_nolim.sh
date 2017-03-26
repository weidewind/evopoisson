#!/bin/bash
foldername="nolimreal"
perl poissonnolim.pl -p h1 --state syn --simnumber 10000 --no_neighbour_changing --mutnum_control 0  --maxmem 120000000 --output $foldername
perl poissonnolim.pl -p h3 --state nsyn --simnumber 10000 --mutnum_control 0  --maxmem 120000000 --output $foldername
perl poissonnolim.pl -p h3 --state syn --simnumber 10000 --no_neighbour_changing --mutnum_control 0  --maxmem 120000000 --output $foldername
perl single_sites_nolim.pl -p h1 --state syn  --mutnum_control 0 --no_neighbour_changing --output $foldername &
perl single_sites_nolim.pl -p h3 --state nsyn  --mutnum_control 0  --output $foldername &
perl single_sites_nolim.pl -p h3 --state syn  --mutnum_control 0 --no_neighbour_changing  --output $foldername &
wait
perl poissonnolim.pl -p n1 --state nsyn --simnumber 10000 --mutnum_control 0  --maxmem 120000000 --output $foldername
perl poissonnolim.pl -p n1 --state syn --simnumber 10000 --no_neighbour_changing --mutnum_control 0  --maxmem 120000000 --output $foldername
perl poissonnolim.pl -p n2 --state nsyn --simnumber 10000 --mutnum_control 0  --maxmem 120000000 --output $foldername
perl poissonnolim.pl -p n2 --state syn --simnumber 10000 --no_neighbour_changing --mutnum_control 0  --maxmem 120000000 --output $foldername
perl single_sites_nolim.pl -p n1 --state nsyn  --mutnum_control 0  --output $foldername &
perl single_sites_nolim.pl -p n1 --state syn  --mutnum_control 0 --no_neighbour_changing  --output $foldername &
perl single_sites_nolim.pl -p n2 --state nsyn  --mutnum_control 0  --output $foldername &
perl single_sites_nolim.pl -p n2 --state syn  --mutnum_control 0 --no_neighbour_changing  --output $foldername &
wait