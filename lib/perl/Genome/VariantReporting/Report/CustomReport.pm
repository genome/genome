package Genome::VariantReporting::Report::CustomReport;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Report::CustomReport {
    is => 'Genome::VariantReporting::Report::WithHeader',
    doc => 'A custom report where you can specify any headers.',
    has_input => {
        header_labels => {
            is => 'Text',
            doc => 'A comma separated list of the headers desired in the report.',
        },
    },
};

sub name {
    return 'custom';
}

sub required_interpreters {
    return qw(position);
}


sub headers {
    my $self = shift;
    return split(/,/, $self->header_labels);
}

1;
