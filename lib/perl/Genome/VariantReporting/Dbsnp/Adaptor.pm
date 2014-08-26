package Genome::VariantReporting::Dbsnp::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Dbsnp::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',

    has_planned_output => [
        joinx_version => { is  => 'Version', },
        info_string => { is => 'Text', },
    ],
    has_provided_output => [
        vcf => {
            is => 'Path',
        }
    ],
};

sub name {
    "dbsnp";
}

sub resolve_provided_attributes {
    my $self = shift;
    $self->vcf($self->provider->get_attribute("dbsnp_vcf"));
}

sub __provided_output_errors__ {
    my $self = shift;
    my $needed = Set::Scalar->new($self->provided_output_names);
    $needed->delete("vcf");
    $needed->insert("dbsnp_vcf");
    return $self->_get_missing_errors(@_, $needed);
}

1;
