#!/bin/bash
# ./selva_pipe.sh -p h1 -c config_noscale.txt -d probedir -u

CONFIG=config_noscale.txt
LENGTH=500
while getopts c:p:d:l:u option
do
		case "${option}" in
			c) CONFIG=${OPTARG};;
			p) PROTEIN=${OPTARG};;
			d) DIR=${OPTARG};;
			l) LENGTH=${OPTARG};;
			u) UPDATE_TREE=true;;
		esac
done

case "${PROTEIN}" in 
	h1) LENGTH=565;;
	h3) LENGTH=566;;
	n1) LENGTH=470;;
	n2) LENGTH=469;;
esac	

OUTNEWICK="${PROTEIN}.${LENGTH}.newick"
MYCONFIG="${PROTEIN}.${LENGTH}.myconfig.txt"
DIRPATH="../data/${DIR}"
export OUTNEWICK
export LENGTH
perl distance_converter.pl -i "${PROTEIN}.l.r.newick" --length $LENGTH -o $OUTNEWICK
cp $CONFIG $MYCONFIG
perl -pi -e 's/^(TREE_FILE).*$/$1 $ENV{OUTNEWICK}/g' $MYCONFIG
perl -pi -e 's/^(LENGTH).*$/$1 $ENV{LENGTH}/g' $MYCONFIG
java -jar Selva-0.1-SNAPSHOT-jar-with-dependencies.jar $MYCONFIG

mkdir $DIRPATH
mv fitnesses.merged.fasta $DIRPATH
mv changetimes.merged.fasta $DIRPATH
mv $MYCONFIG $DIRPATH
mv $OUTNEWICK "${DIRPATH}/${PROTEIN}.l.r.newick"
~/local/emboss/bin/backtranseq -cfile dummycode.cut -sequence allnodes.merged.fasta -outfile "${DIRPATH}/${PROTEIN}.all.fa"

if [ $UPDATE_TREE ]
then 
	echo "Branch lengths will be updated\n"
	cd ..
	perl print_nsyn_tree.pl -p $PROTEIN --input $DIR
	cd selva
	mv "${DIRPATH}/${PROTEIN}.l.r.updated.newick" "${DIRPATH}/${PROTEIN}.l.r.newick"
fi
