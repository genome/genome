package Genome::Model::SomaticCapture;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticCapture {
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
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            valid_values => [0,1,2,3,4],
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
