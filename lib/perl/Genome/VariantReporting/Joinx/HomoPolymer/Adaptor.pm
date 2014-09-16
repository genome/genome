package Genome::VariantReporting::Joinx::HomoPolymer::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

our $MINIMUM_JOINX_VERSION = 1.9;

class Genome::VariantReporting::Joinx::HomoPolymer::Adaptor {
    is => "Genome::VariantReporting::Framework::Component::Adaptor",

    has_planned_output => [
        joinx_version => {
            is => 'Version',
        },
        max_length => {
            is => 'Integer',
        },
        info_string => { 
            is => 'Text', 
        },
    ],
    has_provided_output => [
        homopolymer_list_id => {
            is => 'String',
        },
    ],
};

sub name {
    'homo-polymer';
}


sub __planned_output_errors__ {
    my ($self, $params) = @_;
    my $version = $params->{joinx_version};
    return $self->SUPER::__planned_output_errors__($params), $self->_get_joinx_version_error($version);
}


sub _get_joinx_version_error {
    my ($self, $version) = @_;
    my @errors;

    if ($version < $MINIMUM_JOINX_VERSION) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['joinx_version'],
            desc => "Provided parameter joinx_version: $version not supporting vcf-annotate-homopolymers, use $MINIMUM_JOINX_VERSION and above",
        );
    }

    return @errors;
}

1;
