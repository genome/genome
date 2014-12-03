package Genome::Model::SomaticValidation::Command::CoverageStatsSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::CoverageStatsSummary {
    is => 'Command::V2',
    doc => 'Generate a spreadsheet, tsv file, of coverage metrics for Somatic Validation models.  Duplicate normal samples will only be reported once.',
    has => [
        output_tsv_file => {
            is => 'String',
            doc => 'The output tsv file path to report coverage metrics to',
        },
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'The Somatic Validation models or an expression to resolve the Somatic Validation models.',
        },
    ],
    has_optional => [
        _writer => {
            is => 'Genome::Utility::IO::SeparatedValueWriter',
        }
    ],
};

sub help_detail {
    return "Summarize the coverage stats for all listed somatic validation models.  One line will be output in the tsv file for each sample, region of interest, wingspan and minimum depth filter."
}

sub execute {
    my $self = shift;

    $self->_load_writer;

    my %sample_roi_to_pp;
    for my $model ($self->models) {
        my $roi_name = $model->region_of_interest_set_name;
        my $build = $model->last_succeeded_build;
        unless ($build) {
            die('Failed to find last succeeded build for model: '. $model->id);
        }
        my %sample_name_to_coverage_stats_method = (
            $build->tumor_sample->name => 'coverage_stats_result',
        );
        if ($build->normal_sample) {
            $sample_name_to_coverage_stats_method{$build->normal_sample->name} = 'control_coverage_stats_result';
        }
        for my $sample_name (keys %sample_name_to_coverage_stats_method) {
            my $coverage_stats_result_method = $sample_name_to_coverage_stats_method{$sample_name};
            if ($sample_roi_to_pp{$sample_name}{$roi_name}) {
                # Previously seen sample, either skip or report error when multiple PP used
                if ($sample_roi_to_pp{$sample_name}{$roi_name} eq $build->processing_profile->id) {
                    $self->status_message('Duplicate sample \''. $sample_name .'\' and region of interest \''. $roi_name .'\' allowed, but only reported once in output.');
                } else {
                    $self->error_message('Duplicate sample \''. $sample_name .'\' and region of interest \''. $roi_name .'\' found with different processing profile IDs!.');
                    die($self->error_message);
                }
            } else {
                # First time we have seen this sample and ROI combination
                $sample_roi_to_pp{$sample_name}{$roi_name} = $build->processing_profile->id;
                my $coverage_stats = $build->$coverage_stats_result_method;
                unless ($self->write_sample_coverage($sample_name,$coverage_stats,$roi_name)) {
                    die('Failed to write coverage stats for sample: '. $sample_name);
                }
            }
        }
    }
    $self->_writer->output->close;
    return 1;
}


sub _load_writer {
    my $self = shift;

    my $stats_summary_headers = Genome::Model::Tools::BioSamtools::StatsSummary->default_headers;
    my @headers = ('subject_name','region_of_interest_set_name','wingspan',@{$stats_summary_headers} );

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->output_tsv_file,
        separator => "\t",
        headers => \@headers,
    );
    unless ($writer) {
        die('Failed to open output tsv file path: '. $self->output_tsv_file);
    }
    $self->_writer($writer);
    return 1;
}

sub write_sample_coverage {
    my $self = shift;
    my ($sample_name,$coverage_stats,$roi_name) = @_;

    my $writer = $self->_writer;

    my $coverage_stats_summary_hash_ref = $coverage_stats->coverage_stats_summary_hash_ref();
    my @wingspans = sort {$a <=> $b} keys %{$coverage_stats_summary_hash_ref};
    for my $wingspan (@wingspans) {
        my $wingspan_coverage_stats_summary_hash_ref = $coverage_stats_summary_hash_ref->{$wingspan};
        for my $min_depth (sort { $a <=> $b } keys %{$wingspan_coverage_stats_summary_hash_ref}) {
            my $data = $wingspan_coverage_stats_summary_hash_ref->{$min_depth};
            $data->{subject_name} = $sample_name;
            $data->{region_of_interest_set_name} = $roi_name;
            $data->{wingspan} = $wingspan;
            $writer->write_one($data);
        }
    }
    return 1;
}

1;
