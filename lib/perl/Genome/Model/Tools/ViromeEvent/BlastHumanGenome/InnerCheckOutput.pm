package Genome::Model::Tools::ViromeEvent::BlastHumanGenome::InnerCheckOutput;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastHumanGenome::InnerCheckOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
    #doesn't use fasta file or sample file -- refactor
    has => [
        file_to_run => {
	    is => 'String',
	    doc => 'file to check and re-submit if necessary',
	    is_input => 1,
	}
    ],
};

sub help_brief {
    "Tool to run human genomic blast for file";
}

sub help_detail {
    "Tool to run human genomic blast for file";
}

sub execute {
    my $self = shift;

    unless ( $self->run_blast_for_stage( 'hg_blast' ) ) {
	$self->log_event( "Failed to run human genomic blast for file: ".$self->file_to_run );
	return;
    }
    return 1;
}

1;
