#!/usr/bin/perl

use File::Spec;
use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use Mutmapnolim (realdata_exists, check_realdata_restriction);
use Getopt::Long;
use Getopt::ArgvFile;
use File::Path qw(make_path remove_tree);
use IPC::System::Simple qw(capture); 



my $protein = "h3";
my $state = 'nsyn';
my $input = '';
my $output = '';	# option variable with default value
my $tag = '';
my $muts = "278_INTNODE4195,209_INTNODE4241,209_INTNODE4201";
my $jpeg;
GetOptions (	'protein=s' => \$protein,
		'state=s' => \$state,
		'input=s' => \$input,
		'output=s' => \$output,
		'tag=s' => \$tag,
		'muts=s' => \$muts,
		'jpeg' => \$jpeg,
	);

my @muts = split(/,/, $muts);
my $mutmap = Mutmapnolim->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state});
my @tree_files = $mutmap->print_subtree_with_mutations(\@muts, $tag);
if ($jpeg){
	eval {
		foreach my $file(@tree_files){
			my $jpg = $file.".jpg";
			my $command = "java -jar ~/lib/FigTree_v1.4.3/lib/figtree.jar -graphic JPEG $file $jpg";
			my $logs = capture($command);
			print $logs."\n";
		}
	};
	if (my $exception = $@){
		print STDERR "figtree error: $exception\n";
		exit 3;
	}
}
