package Genome::Model::RnaSeq::Command::BamQc;

use strict;
use warnings;

use File::Path qw(rmtree);

use Genome;

class Genome::Model::RnaSeq::Command::BamQc {
    is => ['Command::V2'],
    has_input => [
        build_id => {
            is => 'Number',
            doc => 'ID of the RnaSeq build upon which to run BamQc',
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build::RnaSeq',
            id_by => 'build_id',
            doc => 'The build upon which to run coverage stats',
        },
    ],
    has_param => [
        lsf_queue => {
            default => 'long',
        },
    ],
    has_optional_output => [
        bam_qc_result_id => {
            is => 'Number',
            doc => 'ID of the result from running this step',
        },
    ],
    has_optional => [
        bam_qc_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged::BamQc',
            id_by => 'bam_qc_result_id',
            doc => 'The result from running this step',
        },
    ],
    doc => 'runs BamQc on the bam(s) produced in the alignment step',
};

sub sub_command_category { 'pipeline steps' }

sub shortcut {
    my $self = shift;

    my %params = $self->params_for_result;
    my $result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get_with_lock(%params);

    if($result) {
        $self->status_message('Using existing result ' . $result->__display_name__);
        return $self->link_result_to_build($result);
    } else {
        return;
    }
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    my %params = (
        $self->params_for_result,
        log_directory => $build->log_directory,
    );

    my $result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get_or_create(%params);

    $self->link_result_to_build($result);


    return 1;
}

sub params_for_result {
    my $self = shift;
    my $build = $self->build;
    my $pp = $build->processing_profile;

    my $alignment_result = $build->alignment_result;

    unless ($alignment_result) {
        die $self->error_message('No alignment result found for build: '. $build->id);
    }

    return (
        alignment_result_id => $alignment_result->id,
        picard_version      => $pp->picard_version || Genome::Model::Tools::Picard->default_picard_version,
        samtools_version    => $pp->samtools_version || Genome::Model::Tools::Sam->default_samtools_version,
        fastqc_version      => $pp->fastqc_version || Genome::Model::Tools::Fastqc->default_fastqc_version,
        samstat_version     => $pp->samstat_version || Genome::Model::Tools::SamStat::Base->default_samstat_version,
        error_rate_version  => Genome::Model::Tools::BioSamtools::ErrorRate->default_errorrate_version,
        error_rate          => $pp->calculate_error_rate || 0,
        read_length         => 1,
        test_name           => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    my $link = join('/', $build->data_directory, 'bam-qc');
    my $label = join('_', 'bam-qc');
    Genome::Sys->create_symlink($result->output_dir, $link);
    $result->add_user(label => $label, user => $build);

    $self->bam_qc_result($result);

    return 1;
}

1;
