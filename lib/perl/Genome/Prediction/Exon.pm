package Genome::Prediction::Exon;

use strict;
use warnings;

use Genome;
use Carp 'confess';

class Genome::Prediction::Exon {
    type_name => 'exon',
    schema_name => 'files',
    data_source => 'Genome::DataSource::Predictions::Exons',
    id_by => [
        exon_name => { is => 'Text' },
        directory => { is => 'Path' },
    ],
    has => [
        start => { is => 'Number' },
        end => { is => 'Number' },
        strand => { is => 'Text' },
        score => { is => 'Text' },
        five_prime_overhang => { is => 'Number' },
        three_prime_overhang => { is => 'Number' },
        transcript_name => { is => 'Text' },
        gene_name => { is => 'Text' },
        sequence_name => { is => 'Text' },
        sequence_string => { is => 'Text' },
        transcript => {
            calculate_from => ['directory', 'transcript_name'],
            calculate => q|
                my ($transcript) = Genome::Prediction::Transcript->get(directory => $directory, transcript_name => $transcript_name);
                return $transcript;
            |,
        },
        coding_gene => {
            calculate_from => ['directory', 'gene_name'],
            calculate => q|
                my ($gene) = Genome::Prediction::CodingGene->get(directory => $directory, gene_name => $gene_name);
                return $gene
            |,
        },
        length => {
            calculate_from => ['start', 'end'],
            calculate => q{ return abs($start - $end) + 1; },
        },
    ],
};

1;
