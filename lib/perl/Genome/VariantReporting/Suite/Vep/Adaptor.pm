package Genome::VariantReporting::Suite::Vep::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Suite::Vep::Adaptor {
    is => "Genome::VariantReporting::Framework::Component::Adaptor",

    has_planned_output => [
        ensembl_version => {
            is => 'String',
            doc => 'Version of ensembl database to use',
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
        custom_annotation_tags => {
            is => 'String',
            is_many => 1,
            is_optional => 1,
            doc => "Custom annotation tags to be used for the --custom option(s)"
        },
        joinx_version => {
            is => 'String',
            doc => 'The joinx version to use for sorting'
        },
        plugins_version => {
            is => 'String',
            doc => 'Version of the vepplugins package to use',
        },
        feature_list_ids => {
            is => 'HASH',
            doc => 'A hash keyed on INFO TAG with values of FeatureList IDs',
            is_translated => 1,
        },
        reference_fasta => {
            is => 'Path',
            is_translated => 1,
            doc => 'Reference fasta file',
        },
    ],
    doc => 'Annotate vcf with results from Ensembl VEP (Variant Effect Predictor)',
};

sub name {
    'vep';
}

1;
