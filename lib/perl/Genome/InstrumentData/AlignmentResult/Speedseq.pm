package Genome::InstrumentData::AlignmentResult::Speedseq;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::AlignmentResult::Speedseq {
    is => 'Genome::InstrumentData::AlignmentResult',
    has => [
    ],
};

sub _run_aligner {
    my $self = shift;

    #Run get_bam_file to fake revivify the per-lane bams
    $self->_prepare_output_directory;
    $self->get_bam_file;
    return 1;
}

#Use merged bam for header since we don't have an original per-lane bam
sub source_bam_path_for_header {
    my $self = shift;
    return $self->get_merged_bam_to_revivify_per_lane_bam;
}

#Don't create flagstat during revivification - postprocess_bam_file will take
#care of that
sub create_bam_flagstat_and_revivify {
    my ($self, $merged_bam, $revivified_bam) = @_;

    my $cmd = Genome::InstrumentData::AlignmentResult::Command::RecreatePerLaneBam->create(
        merged_bam          => $merged_bam,
        per_lane_bam        => $revivified_bam,
        instrument_data_id  => $self->read_and_platform_group_tag_id,
        samtools_version    => $self->samtools_version,
        picard_version      => $self->picard_version,
        bam_header          => $self->bam_header_path,
        include_qc_failed   => 1,
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
    # Always load from the database, since other merged results may have committed since we updated the UR cache
    my @results = Genome::InstrumentData::AlignmentResult::Merged::Speedseq->load(
        'inputs.value_id' => $self->instrument_data_id,
        test_name => $self->test_name,
    );
    return $self->filter_non_database_objects(@results);
}

sub fillmd_for_sam { return 0; }

sub create_BAM_in_staging_directory {
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

    return $self->SUPER::_check_read_count($filtered_bam_rd_ct);
}

sub _promote_data {
    my $self = shift;

    Genome::Sys->copy_file($self->revivified_alignment_bam_file_path . '.flagstat', File::Spec->join($self->temp_staging_directory, 'all_sequences.bam.flagstat'));
    Genome::Sys->copy_file($self->revivified_alignment_bam_file_path . '.md5', File::Spec->join($self->temp_staging_directory, 'all_sequences.bam.md5'));

    return $self->SUPER::_promote_data;
}

1;
