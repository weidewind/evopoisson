package Textbits;
use File::Path qw(make_path);

	sub concat {
		my $site = shift;
		my $node = shift;
		return $site.":".$node;
	}
	
	sub cleave {
		my $site_node = shift;
		return split(/:/, $site_node);
	}
	
	sub iterationFiles {
		my $dirname = shift;
		make_path ($dirname);
		opendir(DH, $dirname);
		my @files = grep { /.*_[0-9]+$/ }readdir(DH); 
		unless (scalar @files > 0){print "No simulation files found in folder $dirname\n";}
		closedir(DH);
		return @files;
	}