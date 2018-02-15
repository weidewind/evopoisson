import os
import sys
import dendropy
import ete3
from dendropy.simulate import treesim
from ete3 import Tree
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-o", "--output", dest="output",
                  help="output file")
parser.add_option("-n", "--num", dest="number", type="int",
                  help="number of tips")
parser.add_option("-l", "--length",
                  type="int", dest="length",
                  help="number of branches inside one line")

(options, args) = parser.parse_args()

outputfolder = "output/little_ksu/article/scheme/"
outputpath = os.path.join(outputfolder,options.output)
print (outputpath)

newick = "("
for i in range(1, options.number
               ):
    newick = newick+"(:0.1):0.1,"
newick = newick+"(:0.1):0.1):0.1;"

for i in range(1,options.length):
    newick = newick.replace('(:', '((:0.1):')
splitted = str.split(newick, ':')
namednewick = splitted.pop(0)
i = 1
for s in splitted:
   if (i == len(splitted)):
      namednewick = namednewick+"root"+":"+s
   else:
      namednewick = namednewick+str(i)+":"+s
   i += 1

if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
f = open(outputpath, 'w')
f.write(namednewick)
f.close()
