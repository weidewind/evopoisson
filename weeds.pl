#!/usr/bin/perl
use Weeds;
use Getopt::Long;
use Getopt::ArgvFile;
use File::Spec;


my $input = "/export/home/popova/workspace/evopoisson/output/h1test_no_mutnum_branch_skipnoskip_mutscontrolled_fixed/nsyn/maxpath_not_subtracted/h1";
my $fails = 0.5;

GetOptions (	
		'input=s' => \$input,
		'fails=s'=> \$fails,
	);

my $weeds = Weeds->new($input);
$weeds->printWeeds("/export/home/popova/workspace/evopoisson/testfiles/weeds0");
$weeds->worstWeeds({fails_threshold => $fails})->printWeeds("/export/home/popova/workspace/evopoisson/testfiles/weeds");
my $ww = Weeds->readWeeds("/export/home/popova/workspace/evopoisson/testfiles/weeds");
$ww->printWeeds("/export/home/popova/workspace/evopoisson/testfiles/weeds2");