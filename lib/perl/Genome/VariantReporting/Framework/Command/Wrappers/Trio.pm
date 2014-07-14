package Genome::VariantReporting::Framework::Command::Wrappers::Trio;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::Command::Wrappers::Trio {
    is => 'Command::V2',
    has_input => [
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
        },
        snvs_plan_file => {
            is => 'File',
        },
        indels_plan_file => {
            is => 'File',
        },
        output_directory => {
            is => 'Path',
            is_output => 1,
        },
        d0_sample => {
            is => 'Genome::Sample',
            doc => 'Discovery sample',
        },
        d30_sample => {
            is => 'Genome::Sample',
            doc => 'Validation sample',
        },
    ],
};

sub execute {
    my $self = shift;
    my @model_pairs = $self->get_model_pairs;
    for my $model_pair (@model_pairs) {
        Genome::Sys->create_directory($model_pair->output_dir);
        Genome::Sys->create_directory($model_pair->reports_directory);
        Genome::Sys->create_directory($model_pair->logs_directory);
        my $resource_file = $model_pair->generate_resource_file;
        $self->run_reports($model_pair);
    }
    return 1;
}

sub get_model_pairs {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Command::Wrappers::ModelPairFactory->create(
        models => [$self->models],
        d0_sample => $self->d0_sample,
        d30_sample => $self->d30_sample,
        output_dir => $self->output_directory,
    );
    return $factory->get_model_pairs;
}

sub run_reports {
    my ($self, $model_pair) = @_;

    for my $variant_type(qw(snvs indels)) {
        my $plan_accessor = join('_', $variant_type, "plan_file");
        my %params = (
            input_vcf => resolve_input_vcf($model_pair, $variant_type),
            variant_type => $variant_type,
            output_directory => $model_pair->reports_directory,
            plan_file => $self->$plan_accessor,
            resource_file => $model_pair->resource_file,
            log_directory => $model_pair->logs_directory,
        );
        Genome::VariantReporting::Framework::Command::CreateReport->execute(%params);
    }
}

sub resolve_input_vcf {
    my ($model_pair, $variant_type) = @_;
    $model_pair->{discovery}->get_detailed_vcf_result($variant_type)->get_vcf($variant_type);
}

1;

