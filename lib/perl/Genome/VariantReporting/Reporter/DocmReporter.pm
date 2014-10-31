package Genome::VariantReporting::Reporter::DocmReporter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(many_samples_available_fields);

class Genome::VariantReporting::Reporter::DocmReporter {
    is => [ 'Genome::VariantReporting::Reporter::WithHeader', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    has => [
    ],
};

sub name {
    return 'docm';
}

sub requires_interpreters {
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
    push @headers, many_samples_available_fields([$self->sample_names]);

    return @headers;
}

1;
