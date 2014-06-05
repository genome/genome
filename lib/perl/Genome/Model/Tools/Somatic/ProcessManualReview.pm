package Genome::Model::Tools::Somatic::ProcessManualReview;

use strict;
use warnings;

use strict;
use warnings;
use File::Basename;
use File::Spec;
use Genome
  ; # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Somatic::ProcessManualReview {
	is  => 'Command',
	has => [
		variant_file =>
		  { is => 'String', is_optional => 0, doc => "variant file" },
		reviewed_file =>
		  { is => 'String', is_optional => 0, doc => "variant-reviewed file" },
		statuses => {
			is          => 'String',
			is_optional => 1,
			doc =>
"statuses to select from the reviewed file (default is S; multiple can be separated by a comma)"
		},
		output_file => {
			is          => 'String',
			is_optional => 1,
			doc =>
"optional, if not provided the output will be placed in variant file's path with a suffix _reviewed_status1-statusN"
		},
	],
};

sub help_brief {
"Displays variants from the variants file that have been manually inspected/curated in the reviewed file"
	  ,;
}

sub execute {
	my $self    = shift;
	my $dirName = dirname(__FILE__);
	my $command = join " ", "java -jar",
	  File::Spec->catfile( $dirName, "ProcessManualReview.jar" );

	$command .= ' --variant-file ';
	$command .= $self->variant_file;

	$command .= ' --reviewed-file ';
	$command .= $self->reviewed_file;

	if ( defined $self->statuses ) {
		$command .= ' --statuses ';
		$command .= $self->statuses;
	}

	if ( defined $self->output_file ) {
		$command .= ' --output-file ';
		$command .= $self->output_file;
	}

	# print($command);
	#	print("\n");
	system($command);
	1;
}

