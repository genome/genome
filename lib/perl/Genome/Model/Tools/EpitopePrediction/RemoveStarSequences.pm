package Genome::Model::Tools::EpitopePrediction::RemoveStarSequences;

use strict;
use warnings;
use Bio::SeqIO;

class Genome::Model::Tools::EpitopePrediction::RemoveStarSequences {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    doc => "Outputs a FASTA file after removing *s from the input Fasta sequence",
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Input Fasta format file',
        },
        output_file => {
            is => 'Text',
            doc => 'The output FASTA file after removing star sequences',
        },
    ],
};

sub execute {
    my $self = shift;

    my $input_filename = $self->input_file ;
    my $output_filename= $self->output_file;

    my $in  = Bio::SeqIO->new(
        -file   => $input_filename ,
        -format => 'Fasta'
    );

    my $out = Bio::SeqIO->new(
        -file => ">".$output_filename ,
        -format => 'Fasta'
    );

    while ( my $seq = $in->next_seq() ) {
        my $seq_string = $seq->seq;
        unless ( $seq_string =~ /[^A-Z]/ ) {
            $out->write_seq($seq);
        }
    }

    return 1;
}

1;
