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
        tumor_sample => {
            is => 'Genome::Sample',
            doc => 'Main tumor sample used for discovery',
        },
        additional_sample => {
            is => 'Genome::Sample',
            doc => 'Additional sample to report readcounts on at discovery variant positions',
        },
    ],
};

sub execute {
    my $self = shift;
    my @model_pairs = $self->get_model_pairs;
    for my $model_pair (@model_pairs) {
        $self->run_reports($model_pair);
    }
    $self->run_summary_stats;
    return 1;
}

sub get_model_pairs {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Command::Wrappers::ModelPairFactory->create(
        models => [$self->models],
        d0_sample => $self->tumor_sample,
        d30_sample => $self->additional_sample,
        output_dir => $self->output_directory,
    );
    return $factory->get_model_pairs;
}

sub run_reports {
    my ($self, $model_pair) = @_;

    for my $variant_type(qw(snvs indels)) {
        my $plan_accessor = join('_', $variant_type, "plan_file");
        my %params = (
            input_vcf => $model_pair->input_vcf($variant_type),
            variant_type => $variant_type,
            output_directory => $model_pair->reports_directory($variant_type),
            plan_file => $self->$plan_accessor,
            resource_file => $model_pair->resource_file,
            log_directory => $model_pair->logs_directory($variant_type),
        );
        Genome::VariantReporting::Framework::Command::CreateReport->execute(%params);
    }
}

sub run_summary_stats {
    my $self = shift;
    Genome::Model::SomaticValidation::Command::AlignmentStatsSummary->execute(
        output_tsv_file => File::Spec->join($self->output_dir, "alignment_summary.tsv"),
        models => [$self->models],
    );
    Genome::Model::SomaticValidation::Command::CoverageStatsSummary->execute(
        output_tsv_file => File::Spec->join($self->output_dir, "coverage_summary.tsv"),
        models => [$self->models],
    );
}

1;

