package Genome::DataSource::Predictions::Exons;

use strict;
use warnings;
use Genome;

class Genome::DataSource::Predictions::Exons {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton' ],
};

sub column_order {
    return [qw(
        exon_name
        start
        end
        strand
        score
        five_prime_overhang
        three_prime_overhang
        transcript_name
        gene_name
        sequence_name
        sequence_string
    )];
}

sub sort_order {
    return ['exon_name'];
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
    return $directory . "/exons.csv";
}

1;

