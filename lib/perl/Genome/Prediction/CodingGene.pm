package Genome::Prediction::CodingGene;

use strict;
use warnings;

use Genome;
use Carp 'confess';

class Genome::Prediction::CodingGene {
    type_name => 'coding gene',
    schema_name => 'files',
    data_source => 'Genome::DataSource::Predictions::CodingGenes',
    id_by => [
        gene_name => { is => 'Text' },
        directory => { is => 'Path' },
    ],
    has => [
        fragment => { is => 'Boolean' },
        internal_stops => { is => 'Boolean' },
        missing_start => { is => 'Boolean' },
        missing_stop => { is => 'Boolean' },
        source => { is => 'Text' },
        strand => { is => 'Text' },
        sequence_name => { is => 'Text' },
        start => { is => 'Number' },
        end => { is => 'Number' },
        length => { 
            calculate_from => ['start', 'end'],
            calculate => q{ return abs($start - $end) - 1; },
        },
        transcript => { 
            calculate_from => ['directory', 'gene_name'],
            calculate => q|
                my ($transcript) = Genome::Prediction::Transcript->get(directory => $directory, coding_gene_name => $gene_name);
                return $transcript;
            |,
        },
        protein => {
            calculate_from => ['directory', 'gene_name'],
            calculate => q|
                my ($protein) = Genome::Prediction::Protein->get(directory => $directory, gene_name => $gene_name);
                return $protein;
            |,
        },
    ],
};

1;
