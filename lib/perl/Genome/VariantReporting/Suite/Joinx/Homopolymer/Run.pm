package Genome::VariantReporting::Suite::Joinx::Homopolymer::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

our $MINIMUM_JOINX_VERSION = 1.9;

class Genome::VariantReporting::Suite::Joinx::Homopolymer::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_planned_transient_argument => [
        joinx_version => {
            is  => 'Version',
            doc => 'joinx version to use',
        },
        homopolymer_list_id => {
            is  => 'Text',
            is_translated => 1,
            doc => 'Homopolymer bed file feature list id',
        },
        max_length => {
            is  => 'Integer',
            doc => 'maximum indel length to annotate as in the homopolymer, default is 2',
        },
        info_string => {
            is  => 'Text',
            doc => 'name of per-allele info field to store the annotation, default is HOMP_FILTER',
        }
    ],
    doc => 'Annotate vcf with information from one homopolymer bed file',
};

sub name {
    'homopolymer';
}

sub result_class {
    'Genome::VariantReporting::Suite::Joinx::Homopolymer::RunResult';
}

sub __planned_errors__ {
    my ($self, $params) = @_;
    my $version = $params->{joinx_version};
    return $self->SUPER::__planned_errors__($params), $self->_get_joinx_version_error($version);
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
