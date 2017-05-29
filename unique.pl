#!/usr/bin/perl
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $input;
my $state = "nsyn";
my $no_neighbour_changing;
my $prefix = 200;
GetOptions (	
		'input=s' =>\$input,
		'state=s' =>\$state,
		'no_neighbour_changing' =>\$no_neighbour_changing,
		'prefix=s' =>\$prefix,
	);

my $dirname = File::Spec->catdir(getcwd(), $input); 

opendir(DH, $dirname);
my @dirs = readdir(DH);
closedir(DH);
if (scalar @dirs == 0){
	die "No files or dirs found in $input!\n";
}

my %size_to_dir;

foreach my $di(sort @dirs){
	my $dipath =  File::Spec->catdir($dirname,$di);
	if (-d $dipath && $di =~ /([0-9]+)_fake/){
		my $fakeno = $1;
		
		$dipath =  File::Spec->catdir($dipath,$state,"maxpath_not_subtracted");
		if ($no_neighbour_changing){
			$dipath =  File::Spec->catdir($dipath,"no_neighbour_changing");
		}
		opendir(DH, $dipath);
		print "fake dir found: $dipath \n";
		my @files = readdir(DH);
		closedir(DH);
		if (scalar @files == 0){
			warn "No files found in $dipath!\n";
			next;
		}
		
		
		my ($realdata) = grep(/(.*)_realdata/, @files);	
		next unless defined $realdata;
		my $realdatapath =  File::Spec->catfile($dipath,$realdata);
		my $size = -s $realdatapath;
		$size_to_dir{$size} = $di;
	}
}	


foreach my $size (keys %size_to_dir){
	my $di = $size_to_dir{$size};
	my $newdi = $prefix.$di;
	print "renaming $di to $newdi in $dirname \n";
	my $dipath =  File::Spec->catdir($dirname,$di);
	rename(File::Spec->catdir($dirname,$di), File::Spec->catdir($dirname,$newdi)) || die ( "Error in renaming: ".$! );
}
