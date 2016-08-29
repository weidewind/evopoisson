use Cwd qw(abs_path cwd getcwd);
use lib getcwd(); #adds working directory to @INC
use MutMap;
use Getopt::Long;


	my $protein;
	my $state = 'nsyn';
	my $input = '';
	my $output = '';	# option variable with default value
	my $verbose;
	GetOptions (	'protein=s' => \$protein,
			'state=s' => \$state,
			'input=s' => \$input,
			'output=s' => \$output,
			'verbose'  => \$verbose,
		);

	my $mutmap = MutMap->new({bigdatatag => $input, bigtag => $output, protein => $protein, state => $state});
	return $mutmap;
