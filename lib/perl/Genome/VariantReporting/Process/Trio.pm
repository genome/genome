package Genome::VariantReporting::Process::Trio;

use strict;
use warnings FATAL => 'all';
use Genome;
use JSON qw(from_json);

class Genome::VariantReporting::Process::Trio {
    is => 'Genome::Process',
    has_input => [
        builds => {
            is => 'Genome::Model::Build::SomaticValidation',
            is_many => 1,
        },
        coverage_builds => {
            is => 'Genome::Model::Build::SomaticValidation',
            is_many => 1,
            is_optional => 1,
        },
        tumor_sample => {
            is => 'Genome::Sample',
        },
        followup_sample => {
            is => 'Genome::Sample',
        },
        normal_sample => {
            is => 'Genome::Sample',
        },
    ],
};

sub get_reports_structure {
    my $self = shift;

    my @report_users = $self->result_users('label like' => 'report:%');

    my $reports_structure = {};
    for my $user (@report_users) {
        my $m;
        if ($user->label =~ /report:(.*)/) {
            my $metadata_json = $1;
            $m = from_json($metadata_json);
        }
        $reports_structure->{$m->{category}}->{$m->{roi_name}}->
            {$m->{variant_type}}->{$m->{report_name}} =
            $user->software_result;
    }

    return $reports_structure;
}

sub get_igv_session_results {
    my $self = shift;

    my @report_users = $self->result_users('label like' => 'igv_session:%');

    my $results = {};
    for my $user (@report_users) {
        my $m;
        if ($user->label =~ /igv_session:(.*)/) {
            my $metadata_json = $1;
            $m = from_json($metadata_json);
        }
        $results->{$m->{roi_name}} = $user->software_result;
    }

    return $results;
}

sub symlink_results {
    my $self = shift;
    my $destination = shift;

    my %structure = %{$self->get_reports_structure};
    while (my ($category, $roi_structure) = each %structure) {
        my %roi_structure = %{$roi_structure};
        while (my ($roi_name, $variant_type_structure) = each %roi_structure) {
            my %variant_type_structure = %{$variant_type_structure};
            while (my ($variant_type, $report_name_structure) = each %variant_type_structure) {
                my %report_name_structure = %{$report_name_structure};
                while (my ($report_name, $report) = each %report_name_structure) {
                    my $link_dir = Genome::Sys->create_directory(File::Spec->join($destination,
                            $category, $roi_name, $variant_type, $report_name));
                    Genome::Sys->symlink_directory($report->output_dir, $link_dir)
                }
            }
        }
    }

    my $alignment_stats_summary_result = $self->results(
        'subclass_name' => 'Genome::Model::SomaticValidation::Command::AlignmentStatsSummary');
    my $alignment_stats_dir = Genome::Sys->create_directory(
        File::Spec->join($destination, 'alignment_stats_summary'));
    Genome::Sys->symlink_directory($alignment_stats_summary_result->output_dir,
        $alignment_stats_dir);

    my $coverage_stats_summary_result = $self->results(
        'subclass_name' => 'Genome::Model::SomaticValidation::Command::CoverageStatsSummary');
    my $coverage_stats_dir = Genome::Sys->create_directory(
        File::Spec->join($destination, 'coverage_stats_summary'));
    Genome::Sys->symlink_directory($coverage_stats_summary_result->output_dir,
        $coverage_stats_dir);

    my %igv_session_results = %{$self->get_igv_session_results};
    while (my($roi_name, $result) = each %igv_session_results) {
        my $link_dir = Genome::Sys->create_directory(File::Spec->join($destination,
                'igv_sessions', $roi_name));
        Genome::Sys->symlink_directory($result->output_dir, $link_dir)
    }

    return 1;
}

sub is_cle_verified {
    my $self = shift;

    return 0 unless $self->SUPER::is_cle_verified(@_);

    for my $result ($self->get_cle_input_results) {
        unless($self->result_is_on_cle_disk_group($result)) {
            return 0;
        }
    }
    return 1;
}

sub get_cle_input_results {
    my $self = shift;

    my @results;
    for my $build ($self->builds, $self->coverage_builds) {
        push @results, $build->get_detailed_vcf_result('snvs');
        push @results, $build->get_detailed_vcf_result('indels');
        push @results, $build->merged_alignment_result;
        push @results, $build->control_merged_alignment_result;
    }
    return @results;
}

1;
