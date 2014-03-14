package Genome::Model::Event::Build::ReferenceAlignment::BamQc;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::ReferenceAlignment::BamQc {
    is  => ['Genome::Model::Event'],
    has_transient_optional => [
        _alignment_result => {is => 'Genome::SoftwareResult'}
    ],
    doc => 'runs BamQc on the bam(s) produced in the alignment step',
};

sub lsf_queue {
    return $ENV{GENOME_LSF_QUEUE_BUILD_WORKER};
}

sub bsub_rusage {
    return '';
}

sub shortcut {
    my $self = shift;

    my %params = $self->params_for_result;
    my $result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get_with_lock(%params);

    if ($result) {
        $self->status_message('Using existing result ' . $result->__display_name__);
        return $self->link_result($result);  #add user and symlink dir to alignment result output dir
    } 
    else {
        return;
    }
}

sub execute {
    my $self  = shift;
    my $build = $self->build;
    my $pp    = $build->processing_profile;

    if ($pp->read_aligner_name eq 'imported') {
        $self->warning_message('Skip BamQc step for imported alignment bam');
        return 1;
    }

    my %params = (
        $self->params_for_result,
        log_directory => $build->log_directory,
    );

    my $result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get_or_create(%params);
    $self->link_result($result);

    return 1;
}

sub params_for_result {
    my $self  = shift;
    my $build = $self->build;
    my $pp    = $build->processing_profile;

    unless ($self->_alignment_result) {
        my $instrument_data_input = $self->instrument_data_input;
        my ($align_result) = $build->alignment_results_for_instrument_data($instrument_data_input->value);

        unless ($align_result) {
            die $self->error_message('No alignment result found for build: '. $build->id);
        }
        $self->_alignment_result($align_result);
    }

    my $picard_version = $pp->picard_version;

    if ($picard_version < 1.40) {
        my $pp_picard_version = $picard_version;
        $picard_version = Genome::Model::Tools::Picard->default_picard_version;
        $self->warning_message('Given picard version: '.$pp_picard_version.' not compatible to CollectMultipleMetrics. Use default: '.$picard_version);
    }

    my $instr_data  = $self->instrument_data;
    #read length takes long time to run and seems not useful for illumina/solexa data
    my $read_length = $instr_data->sequencing_platform =~ /^solexa$/i ? 0 : 1;

    my $error_rate = 1;

    if ($pp->can('read_aligner_name')
        and $pp->can('read_aligner_version')
        and defined $pp->read_aligner_name
        and defined $pp->read_aligner_version
        and ($pp->read_aligner_name eq 'bwamem')
        and ($pp->read_aligner_version =~ /^0\.7\.(5a|7)$/)
    ) {
        $error_rate = 0;
    }

    return (
        alignment_result_id => $self->_alignment_result->id,
        picard_version      => $picard_version,
        samtools_version    => $pp->samtools_version,
        fastqc_version      => Genome::Model::Tools::Fastqc->default_fastqc_version,
        samstat_version     => Genome::Model::Tools::SamStat::Base->default_samstat_version,
        error_rate_version  => Genome::Model::Tools::BioSamtools::ErrorRate->default_errorrate_version,
        error_rate          => $error_rate,
        read_length         => $read_length,
        test_name           => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
}

sub link_result {
    my ($self, $result) = @_;

    my $build = $self->build;
    $result->add_user(label => 'uses', user => $build);

    my $align_result = $self->_alignment_result;
    my $link = join('/', $align_result->output_dir, 'bam-qc-'.$result->id);

    if (-l $link) {
        $self->debug_message("Already found a symlink here.");
    }
    else {
        Genome::Sys->create_symlink($result->output_dir, $link);
        $align_result->_reallocate_disk_allocation;
    }

    my @users = $align_result->users;
    unless (grep{$_->user eq $result}@users) {
        $align_result->add_user(label => 'uses', user => $result);
    }

    return 1;
}

1;
