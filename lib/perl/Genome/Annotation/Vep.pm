package Genome::Annotation::Vep;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Vep {
    is => 'Genome::Annotation::Detail::Command',
    has_input => [
        ensembl_annotation_build_id => {
            is => 'String',
        },
        target_region_set => {
            is => 'Genome::FeatureList',
        },
        segmental_duplications_list => {
            is => 'Genome::FeatureList',
        },
        input_vcf_result => {
            is => 'Genome::SoftwareResult',
        },
        variant_type => { is => 'Text', },
        format => { is => 'String', },
        polyphen => { is => 'String', },
        sift => { is => 'String', },
        condel => { is => 'String', },
        quiet => { is => 'String', },
    ],
    has_optional_output => [
        software_result => {
            is => 'Genome::Annotation::Vep::Result',
            doc => 'The software result created during command execution',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->software_result(Genome::Annotation::Vep::Result->get_or_create($self->input_hash));
    return 1;
}
