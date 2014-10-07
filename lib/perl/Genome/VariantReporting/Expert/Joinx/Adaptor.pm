package Genome::VariantReporting::Expert::Joinx::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Expert::Joinx::Adaptor {
    is => 'Genome::VariantReporting::Framework::Component::Adaptor',
    is_abstract => 1,
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

sub vcf_name {
    my $self = shift;
    return join("_", $self->name, "vcf");
}

sub resolve_provided_attributes {
    my $self = shift;
    $self->vcf($self->provider->get_attribute($self->vcf_name));
}

sub __provided_output_errors__ {
    my $self = shift;
    my $needed = Set::Scalar->new($self->provided_output_names);
    $needed->delete("vcf");
    $needed->insert($self->vcf_name);
    return $self->_get_missing_errors(@_, $needed);
}

1;
