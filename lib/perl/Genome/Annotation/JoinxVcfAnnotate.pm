package Genome::Annotation::JoinxVcfAnnotate;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::JoinxVcfAnnotate {
    is => 'Genome::Annotation::Detail::Command',
    has_input => [
        annotation_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
        },
        input_vcf_result => {
            is => 'Genome::SoftwareResult',
        },
        variant_type => { 
            is => 'Text', 
        },
        info_string => {
            is => 'Text',
        },
        joinx_version => {
            is => 'Text',
        },
    ],
    has_optional_output => [
        software_result => {
            is => 'Genome::Annotation::JoinxVcfAnnotate::Result',
            doc => 'The software result created during command execution',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->software_result(Genome::Annotation::JoinxVcfAnnotate::Result->get_or_create($self->input_hash));
    return 1;
}
