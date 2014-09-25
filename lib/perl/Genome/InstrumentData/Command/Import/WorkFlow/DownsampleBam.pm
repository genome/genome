package Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::WorkFlow::Helpers;

class Genome::InstrumentData::Command::Import::WorkFlow::DownsampleBam {
    is => [qw/ Command::V2 Genome::Model::Tools::Picard::WithDownsampleRatio /],
    has_input => {
        bam_path => {
            is => 'Genome::InstrumentData',
            doc => 'Instrument data to use to create file to reimport.',
        },
    },
    has_output => {
        output_bam_path => {
            is => 'Text',
            calculate_from => [qw/ bam_path /],
            calculate => q| return Genome::InstrumentData::Command::Import::WorkFlow::Helpers->insert_extension_into_bam_path($bam_path, 'downsampled') |,
            doc => 'The path of the downsampled bam.',
        },
    },
    has_optional_transient => { _make_downsample_ratio_required => {}, },
};

sub execute {
    my $self = shift;
    $self->debug_message('Downsample bam...');

    my $downsample_ok = $self->_downsample_bam;
    return if not $downsample_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    my $cleanup_ok = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->remove_paths_and_auxiliary_files($self->bam_path);
    return if not $cleanup_ok;

    $self->debug_message('Downsample bam...done');
    return 1;
}

sub _downsample_bam {
    my $self = shift;

    my $bam_path = $self->bam_path;
    $self->debug_message("Bam path: $bam_path");

    my $output_bam_path = $self->output_bam_path;
    $self->debug_message("Downsampled bam path: $output_bam_path");

    my $downsample_ratio = $self->downsample_ratio;
    $self->debug_message("Downsample ratio: $downsample_ratio");

    my $picard_downsample_cmd = Genome::Model::Tools::Picard::Downsample->create(
        input_file => $bam_path,
        output_file => $output_bam_path,
        downsample_ratio => $downsample_ratio,
        random_seed => 1, # makes bam reproducible
    );
    if ( not $picard_downsample_cmd ) {
        $self->error_message('Failed to create GMT Picard downsample command!');
        return;
    }

    if ( eval{ not $picard_downsample_cmd->execute }) {
        $self->error_message('Failed to execute GMT Picard downsample command!');
        return;
    }

    return 1;
}

sub _verify_read_count {
    my $self = shift;
    $self->debug_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $flagstat = $helpers->load_or_run_flagstat($self->bam_path);
    return if not $flagstat;

    my $downsampled_flagstat = $helpers->load_or_run_flagstat($self->output_bam_path);
    return if not $downsampled_flagstat;

    $self->debug_message('Undownsampled bam read count: '.$flagstat->{total_reads});
    $self->debug_message('Down sampled bam read count: '.$downsampled_flagstat->{total_reads});
    my $expected_read_count = sprintf('%d', $flagstat->{total_reads} * $self->downsample_ratio);
    my $allowed_deviation_pct = .1;
    my $allowed_deviation = sprintf('%d', $expected_read_count * $allowed_deviation_pct);
    my $allowed_min_read_count = $expected_read_count - $allowed_deviation;
    my $allowed_max_read_count = $expected_read_count + $allowed_deviation;
    $self->debug_message('Expected bam read count: '.$expected_read_count);
    $self->debug_message('Allowed deviation: %d%%', $allowed_deviation * 100);
    $self->debug_message('Allowed minimum read count: '.$allowed_min_read_count);
    $self->debug_message('Allowed maximum read count: '.$allowed_max_read_count);

    if ( $downsampled_flagstat->{total_reads} < $allowed_min_read_count 
            or $downsampled_flagstat->{total_reads} > $allowed_max_read_count 
    ) {
        $self->error_message('Downsampled bam read count falls outside the allowed range!');
        return;
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

1;

