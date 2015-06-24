package Genome::InstrumentData::AlignmentResult::Speedseq;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::AlignmentResult::Speedseq {
    is => ['Genome::InstrumentData::AlignmentResult', 'Genome::InstrumentData::AlignmentResult::ReliesOnBwa'],
    has_metric => [
        merged_alignment_result_id => {
            is => 'Text',
        }
    ],
};

sub post_create {
    return 1;
}

sub get_bam_file {
    my $self = shift;

    # Create doesn't do all of the necessary post-processing. This is being
    # delayed until the first time that the bam file is revivified.
    # If we don't have an allocation then this is the first revivification and
    # we will need to do post-processing.
    my $guard = $self->get_bam_lock(__PACKAGE__)->unlock_guard();
    unless ($self->disk_allocations) {
        return $self->_inititalize_revivified_bam;
    }
    else {
        return $self->SUPER::get_bam_file;
    }
}

sub _inititalize_revivified_bam {
    my $self = shift;

    $self->_create_disk_allocation;

    $self->_prepare_working_and_staging_directories;

    $self->debug_message("Preparing the output directory...");
    $self->debug_message("Staging disk usage is " . $self->_staging_disk_usage . " KB");
    my $output_dir = $self->output_dir || $self->_prepare_output_directory;
    $self->debug_message("Alignment output path is $output_dir");

    $self->create_bam_header;

    my $bam_file = $self->SUPER::get_bam_file;

    $self->postprocess_bam_file;

    $self->_compute_alignment_metrics;

    $self->debug_message("Moving results to network disk...");
    $self->_promote_data;

    $self->_reallocate_disk_allocation;

    $self->status_message("Alignment complete.");

    return $bam_file;
}

sub create_bam_header {
    my $self = shift;

    return $self->bam_header_path if -s $self->bam_header_path;

    my $scratch_sam_header_file = $self->prepare_scratch_sam_header_file;

    Genome::Sys->move_file($scratch_sam_header_file, $self->bam_header_path);
    return $self->bam_header_path;
}

sub scratch_sam_file_path {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, 'all_sequences.bam.header');
}

#Don't create flagstat during revivification - postprocess_bam_file will take
#care of that
sub create_bam_flagstat_and_revivify {
    my ($self, $merged_bam, $revivified_bam, $bam_header_path) = @_;

    my %params = (
        merged_bam          => $merged_bam,
        per_lane_bam        => $revivified_bam,
        instrument_data_id  => $self->read_and_platform_group_tag_id,
        samtools_version    => $self->samtools_version,
        picard_version      => $self->picard_version,
        bam_header          => $bam_header_path,
        include_qc_failed   => 1,
    );

    if (-e $self->bam_flagstat_path) {
        $params{comparison_flagstat} = $self->bam_flagstat_path;
    }

    my $cmd = Genome::InstrumentData::AlignmentResult::Command::RecreatePerLaneBam->create(
        %params
    );
    unless ($cmd->execute) {
        die $self->error_message('Failed to execute RecreatePerLaneBam for '.$self->id);
    }
}

sub aligner_params_for_sam_header {
    my $self = shift;
    return "speedseq " . $self->aligner_params;
}

sub get_merged_alignment_results {
    my $self = shift;
    return (Genome::InstrumentData::AlignmentResult::Merged::Speedseq->get($self->merged_alignment_result_id));
}

sub fillmd_for_sam { return 0; }

sub create_BAM_in_staging_directory {
    return 1;
}

sub accepts_bam_input {
    return 1;
}

sub final_staged_bam_path {
    my $self = shift;
    return $self->revivified_alignment_bam_file_path;
}

# Override _check_read_count() from Genome::InstrumentData::AlignmentResult to
# filter reads with secondary alignment flag (0x100) or supplementary
# alignments (0x800)when comparing to the fastq.
sub _check_read_count {
    my ($self, $bam_rd_ct) = @_;

    my $flag = 0x900;
    my $bam_file = $self->final_staged_bam_path;

    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $cmd = "$sam_path view -F $flag -c " . $bam_file;
    my $filtered_bam_rd_ct = `$cmd`;

    $self->debug_message("Overriding _check_read_count: filtering flag $flag from bam read count.");
    $self->debug_message("Actual read count: $bam_rd_ct; filtered read count: $filtered_bam_rd_ct");

    $self->_fastq_read_count($self->determine_input_read_count_from_bam);

    return $self->SUPER::_check_read_count($filtered_bam_rd_ct);
}

sub _create_bam_md5 {
    return 1;
}

sub _create_bam_index {
    return 1;
}

sub _promote_data {
    my $self = shift;

    Genome::Sys->copy_file($self->revivified_alignment_bam_file_path . '.flagstat', File::Spec->join($self->temp_staging_directory, 'all_sequences.bam.flagstat'));

    return $self->SUPER::_promote_data;
}

sub bwa_version {
    my $class = shift;
    my $refindex = shift;

    my $speedseq_version = $refindex->aligner_version;
    my %speedseq_version_to_bwa_version = (
        'test' => '0.7.10',
    );
    my $bwa_version = $speedseq_version_to_bwa_version{$speedseq_version};
    unless ($bwa_version) {
        die $class->error_message('No bwa version assigned to speedseq version (%s)', $speedseq_version);
    }
    return $bwa_version;
}

1;
