#!/usr/bin/perl

use File::Spec;

#my $file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/sites/pandemicH1N1/h1/3lzgsurface.ps";
#my $head=17;
#print_groups($file, $head);

#my $file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/sites/pandemicH1N1/h1/3lzgsurfaceBstrain.ps";
#my $head=344;
#print_groups($file, $head);

my $file = "C:/Users/weidewind/Documents/CMD/Coevolution/Influenza/sites/pandemicH1N1/n1/1-0-1506084778-esp.ps";
my $head=81;
print_groups($file, $head);

sub print_groups {
	my $file = shift;
	my $head = shift;
	open FILE, "<$file" or die "Cannot open file $file: $!\n";
	my @lightblue;
	my @darkblue;
	my $counter=1+$head;
	my $start = $counter;
	
	while (<FILE>){
		if ($_=~ /.*acc\).*/){
			my $str = <FILE>;
			$str = <FILE>;
			chomp $str;
			my $seq = substr($str,index($str,'(')+1);
			for (my $i = 0; $i < length($seq); $i++){
				if (substr($seq, $i, 1) eq 'c'){ push @lightblue, $counter+$i; }
			}
	
			$str = <FILE>;
			$str = <FILE>;
			$str = <FILE>;
			chomp $str;
			$seq = substr($str,index($str,'(')+1);
			for (my $i = 0; $i < length($seq); $i++){
				if (substr($seq, $i, 1) eq 'c'){ push @darkblue, $counter+$i; }
			}
	
			$str = <FILE>;
			$str = <FILE>;
			$str = <FILE>;
			chomp $str;
			$seq = substr($str,index($str,'(')+1);
			$counter += length($seq);
	
			
		}
	}
	close FILE;
	
	my @thicksurf = @{\@darkblue};
	push @thicksurf, @lightblue;
	my %thicksurf = map {$_ => 1} @thicksurf;
	print "thicksurface ";
	foreach my $ts (@thicksurf){
		print $ts.",";
	}
	print "\n";
	my @all = ($start..$counter-1);
	my @internal = grep {not $thicksurf{$_}} @all; 
	print "internal ";
	foreach my $i (@internal){
		print $i.",";
	}
	print "\n";
	
	#foreach my $lb (@lightblue){
	#	print $lb.",";
	#}
	print "\n";
	print "surface ";
	foreach my $db (@darkblue){
		print $db.",";
	}
	print "\n";
}