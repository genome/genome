package Genome::Model::Tools::Joinx::VcfAnnotationSpec;

use Genome;
use strict;
use warnings;

class Genome::Model::Tools::Joinx::VcfAnnotationSpec {
    has => [
        annotation_file => {
            is => 'Text',
            doc => 'Vcf File containing annotation',
        },
        info_fields => {
            is => 'Text',
            doc => 'Field ids to embed from the annotation VCF.',
            is_optional => 1,
        },
    ],
    has_optional => [
        identifiers => {
            is => 'Boolean',
            default => 1,
            doc => 'copy identifiers from the annotation file',
        },
        info => {
            is => 'Boolean',
            default => 1,
            doc => 'copy information from info fields from the annotation file',
        },
    ]
};

1;
