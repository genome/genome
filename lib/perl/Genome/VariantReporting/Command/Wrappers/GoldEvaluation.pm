package Genome::VariantReporting::Command::Wrappers::GoldEvaluation;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;

class Genome::VariantReporting::Command::Wrappers::GoldEvaluation {
    is => 'Command::V2',
    has_input => {
        model => {
            is => 'Genome::Model::SomaticValidation',
        },
        gold_sample_name => {
            is => 'Text',
            doc => 'The name of the sample in the VCF used to retrieve genotypes.',
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
        },
        snvs_input_vcf => {
            is => 'Path',
        },
        indels_input_vcf => {
            is => 'Path',
        },
    },
};

sub execute {
    my $self = shift;

    my $model = $self->model;

    my $model_pair;
    if ($self->is_single_bam($model)) {
        # Germline
        $model_pair = Genome::VariantReporting::Command::Wrappers::SingleModel->create(
            discovery => $model->last_succeeded_build,
            base_output_dir => $self->output_directory,
            plan_file_basename => 'gold_germline_report_TYPE.yaml',
            gold_sample_name => $self->gold_sample_name,
        );
    } else {
        #Somatic
        $model_pair = Genome::VariantReporting::Command::Wrappers::ModelPair->create(
            discovery => $model->last_succeeded_build,
            base_output_dir => $self->output_directory,
            plan_file_basename => 'gold_somatic_report_TYPE.yaml',
            gold_sample_name => $self->gold_sample_name,
        );
    }
    for my $variant_type (qw(snvs indels)) {
        my $vcf_method = $variant_type .'_input_vcf';
        my %params = (
            input_vcf => $self->$vcf_method,
            variant_type => $variant_type,
            output_directory => $model_pair->reports_directory($variant_type),
            plan_file => $model_pair->plan_file($variant_type),
            translations_file => $model_pair->translations_file,
            log_directory => $model_pair->logs_directory($variant_type),
        );
        Genome::VariantReporting::Command::CreateReport->execute(%params);
    }
    return 1;
};

sub is_valid {
    my $self = shift;

    if (my @problems = $self->__errors__) {
        $self->error_message('Germline is invalid!');
        for my $problem (@problems) {
            my @properties = $problem->properties;
            $self->error_message("Property " .
                join(',', map { "'$_'" } @properties) .
                ': ' . $problem->desc);
        }
        return;
    }

    return 1;
}

sub is_single_bam {
    my $self = shift;
    my $model = shift;
    return (!defined($model->normal_sample));
}

1;

