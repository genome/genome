package Genome::VariantReporting::Framework::Test::WithHeaderAndSampleNamesReporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::WithHeaderAndSampleNamesReporter {
    is => 'Genome::VariantReporting::Reporter::WithHeaderAndSampleNames',
};

sub name {
    return '__with_header_and_sample_names__';
}

sub requires_interpreters {
    return qw(__with_many_sample_names__);
}

sub headers {
    my $self = shift;

    $self->create_sample_specific_field_names([$self->fields], [$self->sample_names]);
}

sub fields {
    return qw/info/;
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    $self->_output_fh->print(pp($interpretations) . "\n");
}


1;
