package Genome::Prediction::RNAGene;

use strict;
use warnings;

use Genome;
use Carp 'confess';

class Genome::Prediction::RNAGene {
    type_name => 'rna gene',
    schema_name => 'files',
    data_source => 'Genome::DataSource::Predictions::RNAGenes',
    id_by => [
        gene_name => { is => 'Text' },
        directory => { is => 'DirectoryPath' },
    ],
    has => [
        description => { is => 'Text' },
        start => { is => 'Number' },
        end => { is => 'Number' },
        strand => { is => 'Text' },
        source => { is => 'Text' },
        score => { is => 'Number' },
        sequence_name => { is => 'Text' },
        sequence_string => { is => 'Text' },
    ],
    has_optional => [
        accession => { is => 'Text' },
        product => { is => 'Text' },
        codon => { is => 'Text' },
        amino_acid => { is => 'Text' },
    ],
};

1;
