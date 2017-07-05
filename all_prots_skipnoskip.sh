#!/bin/bash

(perl poissonnolim.pl -p h3 --state nsyn --step 0.5 --restrictions 50,100 --onlysim --simnumber 100 --mutnum_control 0 --shuffler_type exp  --output skipnoskip &>output/logs/h3_skipnoskip ;
perl groups_nolim.pl -p h3 --state nsyn --restrictions 50,100 --mutnum_control 0 --output skipnoskip &>output/logs/h3_skipnoskip_groups &
perl single_sites_nolim.pl -p h3 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/h3_skipnoskipsites) &
(perl poissonnolim.pl -p h1 --state nsyn --step 0.5 --restrictions 50,100 --onlysim --simnumber 100 --mutnum_control 0 --shuffler_type exp  --output skipnoskip &>output/logs/h1_skipnoskip ;
perl groups_nolim.pl -p h1 --state nsyn --restrictions 50,100 --mutnum_control 0 --output skipnoskip &>output/logs/h1_skipnoskip_groups &
perl single_sites_nolim.pl -p h1 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/h1_skipnoskipsites) &
(perl poissonnolim.pl -p n2 --state nsyn --step 0.5 --restrictions 50,100 --onlysim --simnumber 100 --mutnum_control 0 --shuffler_type exp  --output skipnoskip &>output/logs/n2_skipnoskip ;
perl groups_nolim.pl -p n2 --state nsyn --restrictions 50,100 --mutnum_control 0 --output skipnoskip &>output/logs/n2_skipnoskip_groups &
perl single_sites_nolim.pl -p n2 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/n2_skipnoskipsites) &
perl poissonnolim.pl -p n1 --state nsyn --step 0.5 --restrictions 50,100 --onlysim --simnumber 100 --mutnum_control 0 --shuffler_type exp  --output skipnoskip &>output/logs/n1_skipnoskip
perl groups_nolim.pl -p n1 --state nsyn --restrictions 50,100 --mutnum_control 0 --output skipnoskip &>output/logs/n1_skipnoskip_groups &
perl single_sites_nolim.pl -p n1 --state nsyn --step 0.5 --restrictions 50 --mutnum_control 0 --output skipnoskip &>output/logs/n1_skipnoskipsites
wait
