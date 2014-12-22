package Genome::Model::Event::Build::ReferenceAlignment::BamQc;

use strict;
use warnings;

use version 0.77;
use Genome;
use Set::Scalar;

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

        my %segment_info;
        if(defined $self->instrument_data_segment_id) {
            for my $key (qw(instrument_data_segment_id instrument_data_segment_type)) {
                $segment_info{$key} = $self->$key;
            }
        }

        my ($align_result) = $pp->results_for_instrument_data_input($instrument_data_input, %segment_info);

        unless ($align_result) {
            die $self->error_message('No alignment result found for build: '. $build->id);
        }
        $self->_alignment_result($align_result);
    }

    my $picard_version = $self->_select_picard_version($pp->picard_version);


    my $instr_data  = $self->instrument_data;
    #read length takes long time to run and seems not useful for illumina/solexa data
    my $read_length = $instr_data->sequencing_platform =~ /^solexa$/i ? 0 : 1;

    my $error_rate_version = $self->_select_error_rate_version_for_pp($pp);
    my $should_run_error_rate = $self->_should_run_error_rate_for_pp($pp);

    return (
        alignment_result_id => $self->_alignment_result->id,
        picard_version      => $picard_version,
        samtools_version    => $pp->samtools_version,
        fastqc_version      => Genome::Model::Tools::Fastqc->default_fastqc_version,
        samstat_version     => Genome::Model::Tools::SamStat::Base->default_samstat_version,
        error_rate_version  => $error_rate_version,
        error_rate          => $should_run_error_rate,
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

sub _select_picard_version {
    my ($self, $picard_version) = @_;

    my $selected_picard_version = $picard_version;
    #picard versions are a little odd. 1.40 is less than 1.113 but then simple arithmetic doesn't work
    if (version->parse("v$picard_version") < version->parse("v1.40")) {
        $selected_picard_version = Genome::Model::Tools::Picard->default_picard_version;
        $self->warning_message('Requested picard version: '.$picard_version.' is not compatible with CollectMultipleMetrics. Using default: '.$selected_picard_version);
    }

    return $selected_picard_version;
}

sub _select_error_rate_version_for_pp {
    my ($self, $pp) = @_;

    my $error_rate_version = Genome::Model::Tools::BioSamtools::ErrorRate->default_errorrate_version;

    if ($pp->can('read_aligner_name')
        and $pp->can('read_aligner_version')
        and defined $pp->read_aligner_name
        and defined $pp->read_aligner_version
        and ($pp->read_aligner_name =~ /^bwamem/)
    ) {
        my $mem_version = $self->_bwa_mem_version_object($pp->read_aligner_version);
        if($mem_version > $self->_bwa_mem_version_object("0.7.5")) {
            $error_rate_version = '1.0a2';
        }
    }

    return $error_rate_version;
}

sub _should_run_error_rate_for_pp {
    my ($self, $pp) = @_;

    my $aligner_black_listed = ($pp->can('read_aligner_name')
        and $self->_aligner_blacklist->has($pp->read_aligner_name));

    return $aligner_black_listed ? 0 : 1;
}

sub _aligner_blacklist {
    return Set::Scalar->new(qw(bsmap));
}

sub _bwa_mem_version_object {
    my ($self, $mem_version) = @_;
    my ($main_version, $letter_version) = $mem_version =~ /(\d+\.\d+.\d+)([a-z]){0,1}/;
    my $full_version = "v$main_version";
    if(defined $letter_version) {
        $full_version = join("_",$main_version,ord($letter_version) - 96);
    }
    return version->parse($full_version);
}


1;
