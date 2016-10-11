#!/bin/bash
perl print_subtrees.pl --protein h1 --state nsyn -o middle -m 238_INTNODE2422,73_INTNODE2425,151_INTNODE2425,171_INTNODE2429,202_INTNODE2406   --tag likelihood_best_sites_adaptation --jpeg &
perl print_subtrees.pl --protein h1 --state nsyn -o middle -m 205_INTNODE2406,169_INTNODE2434,169_INTNODE2394,144_INTNODE2406,155_INTNODE2428,85_INTNODE2434   --tag likelihood_best_sites_ageing --jpeg &
perl print_subtrees.pl --protein h3 --state nsyn -o middle -m 15_INTNODE4236,153_INTNODE4102 --tag likelihood_best_sites_adaptation --jpeg &
perl print_subtrees.pl --protein h3 --state nsyn -o middle -m 153_INTNODE4232,218_INTNODE3933,149_INTNODE4200,149_INTNODE4230,159_INTNODE4230,19_INTNODE4238,291_INTNODE4258,19_INTNODE4238,69_INTNODE4232,66_INTNODE4232,110_INTNODE4209,175_INTNODE4215 --tag likelihood_best_sites_ageing --jpeg &
perl print_subtrees.pl --protein n1 --state nsyn -o middle -m 382_INTNODE2390,200_INTNODE2384,332_INTNODE2384,336_INTNODE2442 --tag likelihood_best_sites_ageing --jpeg &
perl print_subtrees.pl --protein n2 --state nsyn -o middle -m 466_INTNODE4636,435_INTNODE4571,197_INTNODE4596,248_INTNODE4608,329_INTNODE4636,221_INTNODE4476 --tag likelihood_best_sites_ageing --jpeg &
wait
