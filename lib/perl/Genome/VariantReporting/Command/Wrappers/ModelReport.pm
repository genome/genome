package Genome::VariantReporting::Command::Wrappers::ModelReport;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;

class Genome::VariantReporting::Command::Wrappers::ModelReport {
    is => 'Command::V2',
    has_input => {
        model => {
            is => 'Genome::Model::SomaticValidation',
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
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
            plan_file_basename => 'germline_report_TYPE.yaml',
        );
    } else {
        #Somatic
        $model_pair = Genome::VariantReporting::Command::Wrappers::ModelPair->create(
            discovery => $model->last_succeeded_build,
            base_output_dir => $self->output_directory,
            plan_file_basename => 'somatic_TYPE_report.yaml',
        );
    }
    for my $variant_type(qw(snvs indels)) {
        my %params = (
            input_vcf => $model_pair->input_vcf($variant_type),
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

