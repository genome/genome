package Genome::Model::Tools::ViromeEvent::BlastX_Viral::InnerCheckOutput;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastX_Viral::InnerCheckOutput{
    is => 'Genome::Model::Tools::ViromeEvent',
    #doesn't use fasta file or sample file -- refactor
    has => [
        file_to_run => {
             is => 'String',  
            doc => 'files to rerun repeat masker', 
            is_input => 1,
        }
    ],
};

sub help_brief {
    'Tool to run Viral blastX for an in put file',
}

sub help_detail {
    'Tool to run Viral blastX for an input file',
}

sub execute {
    my $self = shift;

    unless ( $self->run_blast_for_stage( 'blastx_viral' ) ) {
	$self->error_message( "Failed to run Viral blastX for file: ".$self->file_to_run );
	return;
    }
    return 1;
}

1;
