

#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use Mutmapnolim (realdata_exists, check_realdata_restriction, parse_tree);
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);
use IPC::System::Simple qw(capture); 
use Data::Dumper;

my $treefile;
my $outfile = '';	# option variable with default value
my $hazard = 0.2;
my $mutnum = 5;
my $k = 1; # k < 1 - acclimatization, k > 1 - ageing

GetOptions (	#'protein=s' => \$protein,
		'tree=s' => \$treefile,
		'hazard=s' => \$hazard,
		'output=s' => \$outfile,
		'mutnum=i' => \$mutnum,
		'k=s' => \$k,
	);

my $mutnum_control = 2;
my $poisson = 0;
my $constr = Constrains->new(number_of_mutations => $mutnum, stoppers => [], hazard => $hazard);

my $tree = Mutmapnolim::parse_tree($treefile);
	$tree->visit_breadth_first(
		-in   => sub{
			my $node=shift;
			if($node->is_root){
				$node->set_generic('time' => 0);
			}else{
				my $pnode=$node->get_parent;
				my $time=$pnode->get_generic('time');
				$time+=$node->get_branch_length;
				$node->set_generic('time' => $time);
			}
		}
	);	

my @strip_constrs;
push @strip_constrs, $constr;
my @tmp; 
shuffle_muts_on_tree_exp::shuffle_mutations($tree->get_root(),\@strip_constrs,\@tmp, $poisson, $mutnum_control, $k);
print Dumper (\@tmp);
my %events = map {$_ => 1}  @{$tmp[0]};
print Dumper (\%events);
my @subtree_nodes;
print "root node name ".$tree->get_root()->get_name."\n";
#todo: traverse tree and fill in the subtree array 
$tree->visit_depth_first(
		-pre   => sub{
			my $node=shift;
			unless ($node->is_root()){
				if ($node->get_parent->get_generic('prune')){
					$node->set_generic('prune' => 1);
					print("pruned".$node->get_name); 
				}
				push @subtree_nodes, $node unless ($node->get_generic('prune'));
				if(exists $events{$node->get_name}){
					print "yes";
					$node->set_generic('prune' => 1);
					print("pruned".$node->get_name); 
				}
			}
		}
);	

my $obs_hazard;
foreach my $node (@subtree_nodes){
	$obs_hazard += $node->get_branch_length();
}
$obs_hazard = (scalar keys %events)/$obs_hazard;
print "observed hazard $obs_hazard\n";		

my $eventsfilepath = $treefile."_treescheme";
open FILE, ">$eventsfilepath";
print FILE "Ancnode:".$tree->get_root()->get_name()."\n";
print FILE "Events:";
foreach my $node ( @{$tmp[0]}){ #array of nodenames
	print FILE $node.",";
}
print FILE "\n";
print FILE "Subtree:";
foreach my $node (@subtree_nodes){
	print FILE $node->get_name.",";
}
print FILE "\n";
print FILE "observed hazard $obs_hazard (with full event branches, not halves, as in the stats)\n";
print FILE "simulation hazard $hazard\n";
close FILE;



my $command = "xvfb-run python drawTree.py --treefile $treefile --eventfile $eventsfilepath --output $outfile --scale 15 --width 250 --circle_size 2";
my $logs = capture($command);
print $logs."\n";

