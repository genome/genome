package Genome::InstrumentData::AlignmentResult::Speedseq;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::AlignmentResult::Speedseq {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_metric => [
        merged_alignment_result_id => {
            is => 'Text',
        }
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

sub path_for_bam_header_creation {
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

    return $self->SUPER::_promote_data;
}

sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    my $staging_dir = $refindex->temp_staging_directory;

    Genome::Sys->create_symlink($refindex->reference_build->get_sequence_dictionary("sam"), $staging_dir ."/all_sequences.dict" );
    my $bwa_version = $class->bwa_version($refindex->aligner_version);
    my $bwa_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        reference_build_id => $refindex->reference_build_id,
        aligner_name => 'bwa',
        aligner_version => $bwa_version,
        test_name => Genome::Config::get('aligner_index_test_name'),
        users => $refindex->_user_data_for_nested_results,
    );
    for my $filepath (glob($bwa_index->output_dir . "/*")){
        my $filename = File::Basename::fileparse($filepath);
        next if $filename eq 'all_sequences.fa';
        next if $filename eq 'all_sequences.dict';
        Genome::Sys->create_symlink($filepath, $staging_dir . "/$filename");
    }

    $bwa_index->add_user(
        label => 'uses',
        user => $refindex
    );

    return 1;
}


sub bwa_version {
    my $class = shift;
    my $speedseq_version = shift;
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
