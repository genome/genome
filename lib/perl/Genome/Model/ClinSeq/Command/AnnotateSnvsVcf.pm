package Genome::Model::ClinSeq::Command::AnnotateSnvsVcf;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::AnnotateSnvsVcf {
    is => 'Genome::Command::DelegatesToResult',
    has => [
        somatic_build => {
            is => 'Genome::Model::Build::SomaticInterface',
            doc => 'The somatic build with the snvs vcf to annotate',
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
    return 'Genome::Model::ClinSeq::Command::AnnotateSnvsVcf::Result';
}

sub input_hash {
    my $self = shift;

    my %input_hash = $self->SUPER::input_hash;
    $input_hash{input_file}      = $self->somatic_build->snvs_variants_vcf_file;
    $input_hash{annotation_file} = $self->somatic_build->previously_discovered_variations_build->snvs_vcf;

    return %input_hash;
}

1;
