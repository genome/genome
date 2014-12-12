package Genome::VariantReporting::Report::DocmReport;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(many_samples_available_fields);

class Genome::VariantReporting::Report::DocmReport {
    is => [ 'Genome::VariantReporting::Report::WithHeader', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    doc => 'Output readcount information from bam readcount',
};

sub name {
    return 'docm';
}

sub required_interpreters {
    return qw(position many-samples-vaf);
}

sub headers {
    my $self = shift;
    my @headers = qw/
        chromosome_name
        start
        stop
        reference
        variant
    /;
    push @headers, many_samples_available_fields($self);

    return @headers;
}

1;
