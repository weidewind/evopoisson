#!/usr/bin/perl
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;
use Cwd qw(abs_path cwd getcwd);

my $input;
GetOptions (	
		'input=s' =>\$input,
	);

my $dirname = File::Spec->catdir(getcwd(), $input); 
opendir(DH, $dirname);
my @files = readdir(DH);
closedir(DH);
if (scalar @files == 0){
	die "No files found in $input!\n";
}
foreach my $filename(sort @files){
	my $filepath =  File::Spec->catfile($dirname,$filename);
	next if (-d $$filepath);
	next unless ($filename =~ /(.*)_boot_data_([0-9\.]*)_all/);
	my $prot = $1;
	my $outfile = File::Spec->catfile($dirname, $prot."_effects");
	open OUTFILE, ">$outfile" or die "Cannot open $outfile: $!\n";
	print OUTFILE "group,type,stat,bootmean,diff,maxdepth\n";
	close OUTFILE;
}	

foreach my $filename(sort @files){
	my $filepath =  File::Spec->catfile($dirname,$filename);
	next if (-d $$filepath);
	next unless ($filename =~ /(.*)_boot_data_([0-9\.]*)_(.*)/);
	my $prot = $1;
	my $depth = $2;
	my $group = $3;
	open FILE, "<$filepath" or die "Cannot open $filename: $!\n";
	
	my $outfile = File::Spec->catfile($dirname, $prot."_effects");
	open OUTFILE, ">>$outfile" or die "Cannot open $outfile: $!\n";
    my $str = <FILE>;
	my @stattypes = split(/\s+/, $str);
	my $str = <FILE>;
	my @datastat = split(/\s+/, $str);
	while(<FILE>){
		$str = $_;
	}
	my @bootmean = split(/\s+/, $str);
	for (my $i=0; $i < scalar @stattypes; $i++){
		print OUTFILE $group.",".$stattypes[$i].",".$datastat[$i].",".$bootmean[$i].",".($datastat[$i]-$bootmean[$i]).",".$depth."\n";
	}
	close OUTFILE;
	close FILE;
}
