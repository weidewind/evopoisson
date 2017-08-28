#!/bin/bash
called=`basename "$0"`
foldername="output/skipnoskip"
cp $called "${foldername}/${called}"

perl groups_nolim.pl -p h3 --state nsyn --restrictions 50,100 --mutnum_control 0 --output skipnoskip &>output/logs/h3_skipnoskip_groups &
perl single_sites_nolim.pl -p h3 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/h3_skipnoskipsites &
perl single_sites_nolim.pl -p h1 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/h1_skipnoskipsites &
perl single_sites_nolim.pl -p n2 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/n2_skipnoskipsites & 
perl single_sites_nolim.pl -p n1 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/n1_skipnoskipsites &
wait
