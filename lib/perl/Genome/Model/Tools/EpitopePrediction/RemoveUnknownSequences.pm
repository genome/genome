package Genome::Model::Tools::EpitopePrediction::RemoveUnknownSequences;

use strict;
use warnings;
use Bio::SeqIO;
use Genome::Info::CodonToAminoAcid;
use feature "state";

class Genome::Model::Tools::EpitopePrediction::RemoveUnknownSequences {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    doc => "Outputs a FASTA file after removing unknown sequences from the input Fasta sequence",
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Input Fasta format file',
        },
        output_directory => {
            is => 'Text',
            doc => 'Location of the output',
        },
    ],
    has_output => [
        output_file => {
            is => 'Text',
            doc => 'The output FASTA file after removing unknown sequences',
            is_calculated => 1,
            calculate_from => ['output_directory'],
            calculate => q| return File::Spec->join($output_directory, "variant_sequences_filtered.fasta"); |,
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
        if ( is_valid_sequence($seq_string) ) {
            $out->write_seq($seq);
        }
        else {
            $self->warning_message("Sequence for id (%s) contains unknown amino acids: (%s)", $seq->primary_id, $seq_string);
        }
    }

    return 1;
}

sub _generate_valid_sequence_regex {
    my %conversions = %Genome::Info::CodonToAminoAcid::convert;
    delete $conversions{Z};
    my $allowed_amino_acids = join('', keys %conversions);
    return qr/^[$allowed_amino_acids]+$/;
}

sub is_valid_sequence {
    my $seq_string = shift;

    state $valid_sequence_regex = _generate_valid_sequence_regex;
    return $seq_string =~ $valid_sequence_regex;
}

1;
