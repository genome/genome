package Genome::Model::Tools::ViromeEvent::BlastN::InnerCheckOutput;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ViromeEvent::BlastN::InnerCheckOutput{
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
    "Tool to run BlastN for input file";
}
sub help_detail {
    "Tool to run BlastN for input file";
}

sub execute {
    my $self = shift;

    unless ( $self->run_blast_for_stage( 'blast_n' ) ) {
	$self->error_message( "Failed to run BlastN for file: ".$self->file_to_run );
	return;
    }
    return 1;
}

1;
