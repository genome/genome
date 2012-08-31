package Genome::Model::Tools::ViromeEvent::BlastN::OuterCheckOutput;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastN::OuterCheckOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
    has_output=>[
        files_for_blast=> { is => 'ARRAY',  doc => 'array of files for blast n', is_optional => 1},
    ],
};

sub help_brief {
    "Tool that grabs a list of fasta files for running blastN";
}

sub help_detail {
    "Tool that grabs a list of fasta files for running blastN";
}

sub execute {
    my $self = shift;

    unless ( $self->get_files_for_blast( 'blastn' ) ) {
	$self->error_message("Failed to successfully run check blast results");
	return;
    }
    return 1;
}

1;
