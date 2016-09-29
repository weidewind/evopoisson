#!/bin/bash
perl pvalues.pl --protein h1 @cluster.conf &
perl single_sites.pl --protein h1 @cluster.conf &
perl pvalues.pl --protein h3 @cluster.conf &
perl single_sites.pl --protein h3 @cluster.conf &
perl poisson.pl --protein n1 @cluster.conf
perl single_sites.pl --protein n1 @cluster.conf
perl poisson.pl --protein n2 @cluster.conf
perl single_sites.pl --protein n2 @cluster.conf