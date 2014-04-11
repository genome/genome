package Genome::Annotation::Vep;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Vep {
    is => 'Command::V2',
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
            is => 'Genome::Annotation::Readcount::Result',
            doc => 'The software result created during command execution',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->software_result(Genome::Annotation::Vep::Result->get_or_create($self->input_hash));
    return 1;
}

sub input_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(class_name => __PACKAGE__);
    return map {$_->property_name} grep {$_->is_input} @properties;
}

sub input_hash {
    my $self = shift;

    my %hash;
    for my $input_name ($self->input_names) {
        $hash{$input_name} = $self->$input_name;
    }
    return %hash;
}
