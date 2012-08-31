package Genome::Prediction::Transcript;

use strict;
use warnings;

use Genome;
use Carp 'confess';

class Genome::Prediction::Transcript {
    type_name => 'transcript',
    schema_name => 'files',
    data_source => 'Genome::DataSource::Predictions::Transcripts',
    id_by => [
        transcript_name => { is => 'Text' },
        directory => { is => 'Path' },
    ],
    has => [
        coding_gene_name => { is => 'Text' },
        protein_name => { is => 'Text' },
        coding_start => { is => 'Number' },
        coding_end => { is => 'Number' },
        start => { is => 'Number' },
        end => { is => 'Number' },
        strand => { is => 'Text' },
        sequence_name => { is => 'Text' },
        sequence_string => { is => 'Text' },
        coding_gene => {
            calculate_from => ['directory', 'coding_gene_name'],
            calculate => q|
                my ($gene) = Genome::Prediction::CodingGene->get(directory => $directory, gene_name => $coding_gene_name);
                return $gene;
            |,
        },
        protein => {
            calculate_from => ['directory', 'protein_name'],
            calculate => q|
                my ($protein) = Genome::Prediction::Protein->get(directory => $directory, protein_name => $protein_name);
                return $protein;
            |,
        },
        exons => {
            calculate_from => ['directory', 'transcript_name', 'strand'],
            calculate => q|
                my @exons = Genome::Prediction::Exon->get(directory => $directory, transcript_name => $transcript_name);
                @exons = sort { $a->start <=> $b->start } @exons;
                @exons = reverse @exons if $strand eq '-1';
                return @exons;
            |,
        },
    ],
};

# Returns the total length of coding exons
sub spliced_length {
    my $self = shift;
    my $sum = 0;
    for my $exon ($self->exons) {
        $sum += $exon->length;
    }
    return $sum;
}

# Returns the first transcribed exon of the transcript (most five prime)
sub five_prime_exon {
    my $self = shift;
    my @exons = $self->exons;
    return shift @exons;
}
1;
