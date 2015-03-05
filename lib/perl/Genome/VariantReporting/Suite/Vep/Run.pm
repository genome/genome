package Genome::VariantReporting::Suite::Vep::Run;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Vep::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_planned_transient => [
        ensembl_version => {
            is => 'String',
            doc => 'Version of ensembl database to use',
        },
        custom_annotation_tags => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
            doc => "Custom annotation tags to be used for the --custom option(s)",
        },
        feature_list_ids => {
            is => 'HASH',
            is_translated => 1,
            doc => 'A hash keyed on INFO TAG with values of FeatureList IDs',
        },
        reference_fasta => {
            is => 'Path',
            is_translated => 1,
            doc => 'Reference fasta file',
        },
        species => {
            is => 'Text',
            doc => 'Species to use',
        },
        plugins => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
            doc => 'Plugins to use.  These should be formated: PluginName@arg1@arg2,...  If you need to reference the plugin config directory in any of the args, use the placeholder PLUGIN_DIR.  For example, the Condel plugin should be called as follows: Condel@PLUGIN_DIR@b@2',
        },
        joinx_version => {
            is => 'String',
            doc => 'The joinx version to use for sorting',
        },
        plugins_version => {
            is => 'String',
            doc => 'Version of the vepplugins package to use',
        },
    ],
    has_structural_param => [
        lsf_resource => {
            value => q{-R 'select[mem>32000] rusage[mem=32000]' -M 32000000},
        },
    ],
    doc => 'Annotate vcf with results from Ensembl VEP (Variant Effect Predictor)',
};

sub name {
    'vep';
}

sub result_class {
    'Genome::VariantReporting::Suite::Vep::RunResult';
}
