package Genome::Model::Somatic;

use strict;
use warnings;

use Genome;

class Genome::Model::Somatic {
    is  => 'Genome::ModelDeprecated',
    has_param => [
        only_tier_1 => {
            doc => "If set to true, the pipeline will skip ucsc annotation and produce only tier 1 snps",
        },
        min_mapping_quality => {
            doc => "minimum average mapping quality threshold for high confidence call",
        },
        min_somatic_quality => {
            doc => "minimum somatic quality threshold for high confidence call",
        },
        skip_sv => {
            doc => "If set to true, the pipeline will skip structural variation detection",
        },
        sv_detector_version => {
            doc => "Version of the SV detector to use. Reference the Genome::Model::Tools::DetectVariants::SomaticBreakdancer module for currently available versions.",
        },
        sv_detector_params => {
            doc => "Parameters to pass to the SV detector.  For breakdancer, separate params for bam2cfg & BreakDancerMax with a colon. If no parameters are desired, just provide ':'.",
        },
        bam_window_version => {
            doc => "Version to use for bam-window in the copy number variation step.",
        },
        bam_window_params => {
            doc => "Parameters to pass to bam-window in the copy number variation step.",
        },
        sniper_version => {
            doc => "Version to use for bam-somaticsniper for detecting snps and indels.",
        },
        sniper_params => {
            doc => "Parameters to pass to bam-somaticsniper for detecting snps and indels",
        },
        snv_detector_name => {
            doc => "The name of the variant detector to use for snv detection",
        },
        snv_detector_version => {
            doc => "The version of the variant detector to use for snv detection",
        },
        snv_detector_params => {
            doc => "The params to pass to the variant detector to use for snv detection",
        },
        indel_detector_name => {
            doc => "The name of the variant detector to use for indel detection",
        },
        indel_detector_version => {
            doc => "The version of the variant detector to use for indel detection",
        },
        indel_detector_params => {
            doc => "The params to pass to the variant detector to use for indel detection",
        },
        bam_readcount_version => {
            doc => "Version to use for bam-readcount in the high confidence step.",
        },
        bam_readcount_params=> {
            doc => "Parameters to pass to bam-readcount in the high confidence step",
        },
        require_dbsnp_allele_match => {
            doc => "If set to true, the pipeline will require the allele to match during Lookup Variants"  
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3,4 ],#Genome::Model::Tools::Annotate::TranscriptVariants->available_versions ],
        },
    ],
    has_optional => [
        tumor_model_links => {
            is => 'Genome::Model::Link',
            reverse_as => 'to_model',
            where => [ role => 'tumor'],
            is_many => 1,
        },
        tumor_model => {
            is => 'Genome::Model',
            via => 'tumor_model_links',
            to => 'from_model',
        },
        tumor_model_id => {
            is => 'Text',
            via => 'tumor_model',
            to => 'id',
        },
        normal_model_links => {
            is => 'Genome::Model::Link',
            reverse_as => 'to_model',
            where => [ role => 'normal'],
            is_many => 1,
        },
        normal_model => {
            is => 'Genome::Model',
            via => 'normal_model_links',
            to => 'from_model',
        },
        normal_model_id => {
            is => 'Text',
            via => 'normal_model',
            to => 'id',
        },
    ],
};

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

1;
