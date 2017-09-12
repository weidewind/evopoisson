#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); # adds working directory to @INC
use Mutmapnolim;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use Groups;
use Getopt::ArgvFile;
use List::Util;
use POSIX qw(floor ceil);
use List::Util;


my $protein;


GetOptions (	
'protein=s' => \$protein,
);
			
print $protein;