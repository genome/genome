package Genome::VariantReporting::Reporter::BedReporter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Reporter::BedReporter {
    is => 'Genome::VariantReporting::Reporter::WithHeader',
    doc => 'Output variants in bed format',
};

sub can_be_combined {
    return 1;
}

sub combine_parameters {
    return {
        sort_columns => [qw(1 2 3)],
        contains_header => 0,
    };
}

sub name {
    return 'bed';
}

sub required_interpreters {
    return qw(bed-entry);
}

sub headers {
    return qw/
        chromosome_name
        start
        stop
    /;
}

sub print_headers {
    return;
}
1;
