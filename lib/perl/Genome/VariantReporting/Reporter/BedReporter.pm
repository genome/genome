package Genome::VariantReporting::Reporter::BedReporter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Reporter::BedReporter {
    is => 'Genome::VariantReporting::Reporter::WithHeader',
};

sub name {
    return 'bed';
}

sub requires_interpreters {
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
