package Genome::Model::SomaticValidation::Command::CoverageStatsSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::CoverageStatsSummary {
    is => 'Command::V2',
    doc => 'Generate a spreadsheet, tsv file, of coverage metrics for Somatic Validation builds.  Duplicate normal samples will only be reported once.',
    has => [
        output_tsv_file => {
            is => 'String',
            doc => 'The output tsv file path to report coverage metrics to',
        },
        builds => {
            is => 'Genome::Model::Build::SomaticValidation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'The Somatic Validation builds or an expression to resolve the Somatic Validation builds.',
        },
    ],
    has_optional => [
        _writer => {
            is => 'Genome::Utility::IO::SeparatedValueWriter',
        }
    ],
};

sub help_detail {
    return "Summarize the coverage stats for all listed somatic validation builds.  One line will be output in the tsv file for each sample, region of interest, wingspan and minimum depth filter."
}

sub execute {
    my $self = shift;

    $self->_load_writer;

    my %sample_roi_to_pp;
    for my $build ($self->builds) {
        my $roi_name = $build->region_of_interest_set_name;

        my $tumor_sample = $build->tumor_sample;
        if ($sample_roi_to_pp{$tumor_sample->name}{$roi_name}) {
            die('Duplicate builds for tumor sample \''. $tumor_sample->name .'\' region of interest \''. $roi_name .'\' found!.');
        }
        $sample_roi_to_pp{$tumor_sample->name}{$roi_name} = $build->processing_profile->id;
        my $tumor_coverage_stats = $build->coverage_stats_result;
        unless ($self->write_sample_coverage($tumor_sample,$tumor_coverage_stats,$roi_name)) {
            die('Failed to write coverage stats for tumor sample: '. $tumor_sample->name);
        }

        my $normal_sample = $build->normal_sample;
        if ($sample_roi_to_pp{$normal_sample->name}{$roi_name}) {
            if ($sample_roi_to_pp{$normal_sample->name}{$roi_name} eq $build->processing_profile->id) {
                $self->status_message('Duplicate normal sample \''. $normal_sample->name .'\' and region of interest \''. $roi_name .'\' allowed, but only reported once in output.');
            } else {
                die('Duplicate normal sample \''. $normal_sample->name .'\' and region of interest \''. $roi_name .'\' found with different processing profile IDs.  Exiting!');
            }
        } else {
            $sample_roi_to_pp{$normal_sample->name}{$roi_name} = $build->processing_profile->id;
            my $normal_coverage_stats = $build->control_coverage_stats_result;
            unless ($self->write_sample_coverage($normal_sample,$normal_coverage_stats,$roi_name)) {
                die('Failed to write coverage stats for normal sample: '. $normal_sample->name);
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
    my ($sample,$coverage_stats,$roi_name) = @_;

    my $writer = $self->_writer;

    my $coverage_stats_summary_hash_ref = $coverage_stats->coverage_stats_summary_hash_ref();
    my @wingspans = sort {$a <=> $b} keys %{$coverage_stats_summary_hash_ref};
    for my $wingspan (@wingspans) {
        my $wingspan_coverage_stats_summary_hash_ref = $coverage_stats_summary_hash_ref->{$wingspan};
        for my $min_depth (sort { $a <=> $b } keys %{$wingspan_coverage_stats_summary_hash_ref}) {
            my $data = $wingspan_coverage_stats_summary_hash_ref->{$min_depth};
            $data->{subject_name} = $sample->name;
            $data->{region_of_interest_set_name} = $roi_name;
            $data->{wingspan} = $wingspan;
            $writer->write_one($data);
        }
    }
    return 1;
}

1;
