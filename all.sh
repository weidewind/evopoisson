#!/bin/bash
perl poisson.pl --protein h1 --state nsyn --output synresearch --no_leaves @cluster.conf 
perl poisson.pl --protein h3 --state nsyn --output synresearch --no_leaves @cluster.conf 
perl poisson.pl --protein n1 --state nsyn --output synresearch --no_leaves @cluster.conf 
perl poisson.pl --protein n2 --state nsyn --output synresearch --no_leaves @cluster.conf 

#perl single_sites.pl --protein h1 @cluster.conf &
#perl poisson.pl --protein h3 @cluster.conf
#perl single_sites.pl --protein h3 @cluster.conf &
#perl poisson.pl --protein n1 @cluster.conf
#perl single_sites.pl --protein n1 @cluster.conf &
#perl poisson.pl --protein n2 @cluster.conf
#perl single_sites.pl --protein n2 @cluster.conf
