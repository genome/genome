package Genome::Annotation::Vep;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Vep {
    is => 'Genome::Annotation::Detail::Command',
    has_input => [
        ensembl_version => {
            is => 'String',
        },
        feature_list_ids_and_tags => {
            is => 'String',
            is_many => 1,
            doc => 'List of feature lists to be annotated in the 
                    INFO field, along with the tag to be used
                    e.g. 12345:SEGDUP,58676:ROI
                    The id and tag should be separated by a colon',
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
