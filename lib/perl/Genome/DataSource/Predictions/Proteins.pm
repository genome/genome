package Genome::DataSource::Predictions::Proteins;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Predictions::Proteins {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton' ],
};

sub column_order {
    return [qw(
        protein_name
        fragment
        internal_stops
        transcript_name
        gene_name
        sequence_name
        sequence_string
        cellular_localization
        cog_id
        enzymatic_pathway_id
    )];
}

sub sort_order {
    return ['protein_name'];
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
    return $directory . "/proteins.csv";
}

1;

