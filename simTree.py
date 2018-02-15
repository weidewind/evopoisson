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
parser.add_option("-n", "--num", dest="num_tips", type="int",
                  help="number of tips")
parser.add_option("-s", "--star",
                  action="store_true", dest="star", default=False,
                  help="star tree")

(options, args) = parser.parse_args()

outputfolder = "output/little_ksu/article/scheme/"
outputpath = os.path.join(outputfolder,options.output)
print (outputpath)
taxa = dendropy.TaxonNamespace(['Tip']*options.num_tips)
#taxa = dendropy.TaxonNamespace(list(range(1,options.num_tips+1)))
if options.star:
   taxa = dendropy.TaxonNamespace(['']*options.num_tips)
   print(taxa)
   t = treesim.star_tree(taxa)
else:
   t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=options.num_tips,taxon_namespace=taxa)
newick = t.as_string(schema='newick',suppress_rooting=True)
splitted = str.split(newick, ':')
namednewick = splitted.pop(0)
i = 1
for s in splitted:
  # namednewick = namednewick+str(i)+":"+s
   if (i == len(splitted)):
      namednewick = namednewick+"root"+":"+s
   else:
      namednewick = namednewick+str(i)+":"+s
   i += 1
t.print_plot()

if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)
f = open(outputpath, 'w')
f.write(namednewick)
f.close()

sys.stdout = open(outputpath+"plot", 'w')
t.print_plot()



#etetree = Tree(newick)
