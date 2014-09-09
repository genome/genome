package Genome::VariantReporting::Framework::Command::Wrappers::GermlineOnly;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;

class Genome::VariantReporting::Framework::Command::Wrappers::GermlineOnly {
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
    
    my $single_model = Genome::VariantReporting::Framework::Command::Wrappers::SingleModel->create(
        discovery => $model->last_succeeded_build,
        base_output_dir => $self->output_directory,
        plan_file_basename => 'germline_report_TYPE.yaml',
    );
    for my $variant_type(qw(snvs indels)) {
        my %params = (
            input_vcf => $single_model->input_vcf($variant_type),
            variant_type => $variant_type,
            output_directory => $single_model->reports_directory($variant_type),
            plan_file => $single_model->plan_file($variant_type),
            resource_file => $single_model->resource_file,
            log_directory => $single_model->logs_directory($variant_type),
        );
        Genome::VariantReporting::Framework::Command::CreateReport->execute(%params);
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

1;

