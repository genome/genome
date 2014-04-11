package Genome::Annotation::Adaptor;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Adaptor {
    is => 'Command::V2',
    is_abstract => 1,
    has_input => [
        build => {
            is => 'Genome::Model::Build::RunsDV2',
        },
    ],
    has_output => [
        bam_results => {
            is_many => 1,
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        snv_vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
        },
        indel_vcf_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->bam_results($self->resolve_bam_results);
    $self->snv_vcf_result($self->resolve_snv_vcf_result);
    $self->indel_vcf_result($self->resolve_indel_vcf_result);
    $self->annotation_build($self->resolve_annotation_build);
    return 1;
}

sub resolve_snv_vcf_result {
    my $self = shift;
    my $result = eval {$self->build->get_detailed_snvs_vcf_result};
    my $error = $@;
    if ($error) {
        $self->debug_message("No snv result found on build %s:\n%s", $self->build->id, $error);
    }
    return $result;
}

sub resolve_indel_vcf_result {
    my $self = shift;
    my $result = eval {$self->build->get_detailed_indels_vcf_result};
    my $error = $@;
    if ($error) {
        $self->debug_message("No indel result found on build %s:\n%s", $self->build->id, $error);
    }
    return $result;
}
1;

