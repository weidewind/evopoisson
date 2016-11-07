#PBS -S /bin/bash
#PBS -N syn_no_neighbour_changing
#PBS -l mem=30gb
#PBS -l nodes=1:ppn=15
#PBS -m ae
#PBS -M weidewind@gmail.com
#PBS -d /export/home/popova/workspace/evopoisson
#PBS -e /export/home/popova/workspace/evopoisson/logs/syn_no_leaves.err
#PBS -o /export/home/popova/workspace/evopoisson/logs/syn_no_leaves.out

echo Working directory is $PBS_O_WORKDIR
export PERL5LIB=/export/home/popova/perl5/lib/perl5/x86_64-linux:PERL5LIB 
export PERL5LIB=/export/home/popova/perl5/lib/perl5:PERL5LIB
export PERL5LIB=/export/home/popova/.perl/lib/perl5/5.22.1/x86_64-linux:PERL5LIB 
export PERL5LIB=/export/home/popova/.perl/lib/perl5/5.22.1:PERL5LIB
export PERL5LIB=/export/home/popova/.perl/lib/perl5/x86_64-linux:PERL5LIB 
export PERL5LIB=/export/home/popova/.perl/lib/perl5:PERL5LIB
export PERL5LIB=/export/home/popova/perl5/lib/perl5/x86_64-linux:PERL5LIB 
export PERL5LIB=/export/home/popova/perl5/lib/perl5:PERL5LIB 
export PERL5LIB=/export/home/popova/perl5/lib/perl5/x86_64-linux:PERL5LIB 
export PERL5LIB=/export/home/popova/perl5/lib/perl5:PERL5LIB 
~/perl5/perlbrew/perls/perl-5.22.1/bin/perl poisson.pl --protein h1 --state syn --output synresearch --no_neighbour_changing @cluster.conf
echo h1 done 
~/perl5/perlbrew/perls/perl-5.22.1/bin/perl poisson.pl --protein h3 --state syn --output synresearch --no_neighbour_changing @cluster.conf 
echo h3 done 
~/perl5/perlbrew/perls/perl-5.22.1/bin/perl poisson.pl --protein n1 --state syn --output synresearch --no_neighbour_changing @cluster.conf 
echo n1 done 
~/perl5/perlbrew/perls/perl-5.22.1/bin/perl poisson.pl --protein n2 --state syn --output synresearch --no_neighbour_changing @cluster.conf
echo n2 done 
