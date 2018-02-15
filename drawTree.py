from ete3 import Tree, TreeStyle, NodeStyle
from optparse import OptionParser

## example xvfb-run python drawTree.py -t data/little_ksu/h3.l.r.newick -o output/treetestlk.svg -e output/little_ksu/nsyn/skip_stoppers/trees/etetest/h3_treescheme_291_Node_3506
parser = OptionParser()
parser.add_option("-t", "--treefile", dest="treefile",
                  help="newick tree file")
parser.add_option("-e", "--eventfile", dest="eventfile",
                  help="events file")
parser.add_option("-o", "--output", dest="output",
                  help="output file")
parser.add_option("-s", "--scale", dest="scale", type="int",
                  help="8000 for little ksu, 3 for krya")
parser.add_option("-c", "--circle_size", dest="size", type="float",
                  help="20 for little ksu, ? for krya")
parser.add_option("-w", "--width", dest="width", type="int",
                  help="870 for little ksu, ? for krya")
(options, args) = parser.parse_args()
t = Tree(options.treefile, format=1) #flexible with internal node names

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = False
#ts.optimal_scale_level = "full"
ts.scale = options.scale #3 for krya # pixels per branch length unit
ts.branch_vertical_margin = 2
t.ladderize()

efile = open(options.eventfile, "r")
anc = efile.readline().split(':')[1].rstrip('\n')
events = efile.readline().split(':')[1].split(',')
subtree = efile.readline().split(':')[1].split(',')

# Creates an independent node style for each node, which is
# initialized with a red foreground color.
treecolor = "mediumseagreen"
eventscolor = "Black"
for n in t.traverse():
   nstyle = NodeStyle()
   nstyle["fgcolor"] = "red"
   nstyle["size"] = 0

   nstyle["vt_line_width"] = 1
   nstyle["hz_line_width"] = 1
   if n.name in subtree:
      nstyle["hz_line_color"] = treecolor
      nstyle["vt_line_color"] = treecolor
   if n.name in events:   
      nstyle["fgcolor"] = eventscolor
      nstyle["size"] = options.size
      nstyle["vt_line_color"] = eventscolor
   if n.name == anc:
      nstyle["fgcolor"] = treecolor
      nstyle["size"] = options.size
      nstyle["vt_line_color"] = treecolor
   n.set_style(nstyle)

#for e in t.search_nodes(Name="INTNODE1220"):
#	estyle = NodeStyle()
#	estyle["fgcolor"] = "blue"
#	e.set_style(estyle)
	
t.render(options.output + ".pdf", tree_style=ts, dpi=300, w=options.width, units="mm")
t.render(options.output + ".png", tree_style=ts, dpi=300, w=options.width, units="mm")
t.render(options.output + ".svg", tree_style=ts, dpi=300, w=options.width, units="mm")