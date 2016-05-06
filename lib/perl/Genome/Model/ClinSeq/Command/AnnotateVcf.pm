package Genome::Model::ClinSeq::Command::AnnotateVcf;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::AnnotateVcf {
    is => 'Genome::Command::DelegatesToResult',
    has => [
        input_file => {
            is => 'Text',
            doc => 'Vcf File to filter',
        },
        annotation_file => {
            is => 'Text',
            doc => 'Vcf File containing annotation',
        },
    ],
    has_optional_input => [
        info_fields => {
            is => 'Text',
            doc => 'Field ids to embed from the annotation VCF. Use colons to separate multiple field descriptors.',
            #doing the above because UR autosplits on commas with is_many, but joinx uses commas in its field descriptors
        },
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
    ],
};

sub result_class {
    return 'Genome::Model::ClinSeq::Command::AnnotateVcf::Result';
}

1;
