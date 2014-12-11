package Genome::VariantReporting::Report::BedReport;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Report::BedReport {
    is => 'Genome::VariantReporting::Report::WithHeader',
    doc => 'Output variants in bed format',
};

sub can_be_merged {
    return 1;
}

sub merge_parameters {
    return {
        sort_columns => [qw(1 2 3)],
        contains_header => 0,
    };
}

sub file_name {
    return 'report.bed';
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
