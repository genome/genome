package Genome::Model::Tools::ViromeEvent::BlastX_Viral::OuterCheckOutput;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastX_Viral::OuterCheckOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
    has_output=>[
        files_for_blast=> { is => 'ARRAY',  doc => 'array of files for blast n', is_optional => 1},
    ],
};

sub help_brief {
    "Tool to grab list of files for viral blastX";
}

sub help_detail {
    "Tool to grab list of files for viral blastX";
}

sub execute {
    my $self = shift;

    unless ( $self->get_files_for_blast( 'blastx_viral' ) ) {
	$self->error_message("Failed to successfully run check blast results");
	return;
    }
    return 1;
}

1;
