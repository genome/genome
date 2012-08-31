package Genome::DataSource::Predictions::RNAGenes;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Predictions::RNAGenes {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton' ],
};

sub column_order {
    return [qw(
        sequence_name
        gene_name
        description
        start
        end
        strand
        source
        score
        sequence_string
        accession
        product
        codon
        amino_acid
    )];
}

sub sort_order {
    return ['gene_name'];
}

sub delimiter {
    return ",";
}

sub constant_values {
    return ['directory'];
}

sub skip_first_line {
    return 0;
}

sub required_for_get {
    return ['directory'];
}

sub file_resolver {
    my $directory = shift;
    return $directory . "/rna_genes.csv";
}

1;

