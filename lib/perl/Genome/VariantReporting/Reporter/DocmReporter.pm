package Genome::VariantReporting::Reporter::DocmReporter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;

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
    push @headers, $self->_vaf_headers;

    return @headers;
}

sub _vaf_headers {
    my $self = shift;
    Genome::VariantReporting::Suite::BamReadcount::ManySamplesVafInterpreter::available_fields($self);
}

sub vaf_fields {
    my $self = shift;
    Genome::VariantReporting::Suite::BamReadcount::ManySamplesVafInterpreter::vaf_fields($self);
}

1;
