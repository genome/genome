package Genome::DataSource::Predictions::CodingGenes;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Predictions::CodingGenes {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton' ],
};

sub column_order {
    return [qw(
        gene_name
        fragment
        internal_stops
        missing_start
        missing_stop
        source
        strand
        sequence_name
        start
        end
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
    return $directory . "/coding_genes.csv";
}

1;

