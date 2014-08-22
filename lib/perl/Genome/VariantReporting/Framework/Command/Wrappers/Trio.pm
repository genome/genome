package Genome::VariantReporting::Framework::Command::Wrappers::Trio;

use strict;
use warnings;
use Genome;
use File::Basename qw(basename);

class Genome::VariantReporting::Framework::Command::Wrappers::Trio {
    is => 'Command::V2',
    has_input => [
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
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
        normal_sample => {
            is => 'Genome::Sample',
            doc => 'Normal sample',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->run_summary_stats;
    my @model_pairs = $self->get_model_pairs;
    for my $model_pair (@model_pairs) {
        $self->run_reports($model_pair);
    }
    my @roi_directories = map {basename $_} glob(File::Spec->join($self->output_directory, "discovery", "*"));
    for my $roi_directory (@roi_directories) {
        for my $base (Genome::VariantReporting::Framework::Command::Wrappers::ModelPair->report_names) {
            my $discovery_report = File::Spec->join($self->output_directory, "discovery", $roi_directory, $base);
            my $additional_report = File::Spec->join($self->output_directory, "additional", $roi_directory, $base);
            Genome::VariantReporting::PostProcessing::CombineReports->execute(
                reports => [$discovery_report, $additional_report],
                sort_columns => [qw(chromosome_name start stop reference variant)],
                contains_header => 1,
                output_file => File::Spec->join($self->output_directory, "$base-$roi_directory"),
                entry_sources =>  {$discovery_report => 'd0', $additional_report => 'd30'},
            );
        }
    }
    return 1;
}

sub get_model_pairs {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Command::Wrappers::ModelPairFactory->create(
        models => [$self->models],
        d0_sample => $self->tumor_sample,
        d30_sample => $self->additional_sample,
        normal_sample => $self->normal_sample,
        output_dir => $self->output_directory,
    );
    return $factory->get_model_pairs;
}

sub run_reports {
    my ($self, $model_pair) = @_;

    for my $variant_type(qw(snvs indels)) {
        my %params = (
            input_vcf => $model_pair->input_vcf($variant_type),
            variant_type => $variant_type,
            output_directory => $model_pair->reports_directory($variant_type),
            plan_file => $model_pair->plan_file($variant_type),
            resource_file => $model_pair->resource_file,
            log_directory => $model_pair->logs_directory($variant_type),
        );
        Genome::VariantReporting::Framework::Command::CreateReport->execute(%params);
    }
    for my $base ($model_pair->report_names) {
        Genome::VariantReporting::PostProcessing::CombineReports->execute(
            reports => [File::Spec->join($model_pair->reports_directory("snvs"), $base),
                File::Spec->join($model_pair->reports_directory("indels"), $base)],
            sort_columns => [qw(chromosome_name start stop reference variant)],
            contains_header => 1,
            output_file => File::Spec->join($model_pair->output_dir, $base),
            split_indicators => [qw(per_library)],
        );
    }
}

sub run_summary_stats {
    my $self = shift;
    Genome::Model::SomaticValidation::Command::AlignmentStatsSummary->execute(
        output_tsv_file => File::Spec->join($self->output_directory, "alignment_summary.tsv"),
        models => [$self->models],
    );
    Genome::Model::SomaticValidation::Command::CoverageStatsSummary->execute(
        output_tsv_file => File::Spec->join($self->output_directory, "coverage_summary.tsv"),
        models => [$self->models],
    );
}

1;

