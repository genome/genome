package Genome::Model::Tools::EpitopePrediction::FilterSequences;

use strict;
use warnings;
use Bio::SeqIO;
use Genome::Info::CodonToAminoAcid;
use feature "state";

class Genome::Model::Tools::EpitopePrediction::FilterSequences {
    is => 'Genome::Model::Tools::EpitopePrediction::Base',
    doc => "Outputs a FASTA file after removing unknown sequences as well as duplicate sequences from the input FASTA sequence",
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Input FASTA format file',
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

    my %sequences;
    while ( my $seq = $in->next_seq() ) {
        my $seq_string = $seq->seq;
        if (!defined($seq_string)) {
            $self->warning_message("Sequence for id (%s) is undefined", $seq->primary_id);
        }
        elsif ( is_valid_sequence($seq_string) ) {
            if (defined(my $existing_seq_string = $sequences{$seq->primary_id})) {
                if ($existing_seq_string eq $seq_string) {
                    $self->warning_message("Sequence (%s) with id (%s) is a duplicate. Skipping.", $seq_string, $seq->primary_id);
                }
                else {
                    die $self->error_message("Found duplicate entries with id (%s) but different sequences: (%s) and (%s).", $seq->primary_id, $existing_seq_string, $seq_string);
                }
            }
            else  {
                $out->write_seq($seq);
            }
            $sequences{$seq->primary_id} = $seq_string;
        }
        else {
            $self->warning_message("Sequence for id (%s) contains unknown amino acids: (%s). Skipping.", $seq->primary_id, $seq_string);
        }
    }

    return 1;
}

sub _generate_valid_sequence_regex {
    my @valid_amino_acid_codes = Genome::Info::CodonToAminoAcid::valid_amino_acid_codes;
    my $allowed_amino_acids = join('', @valid_amino_acid_codes);
    return qr/^[$allowed_amino_acids]+$/;
}

sub is_valid_sequence {
    my $seq_string = shift;

    state $valid_sequence_regex = _generate_valid_sequence_regex;
    return $seq_string =~ $valid_sequence_regex;
}

1;
