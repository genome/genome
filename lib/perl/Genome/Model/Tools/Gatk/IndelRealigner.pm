package Genome::Model::Tools::Gatk::IndelRealigner;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Gatk::IndelRealigner {
    doc => "Run GATK with the 'IndelRealigner' tool",
    is => 'Genome::Model::Tools::Gatk::Base',
    has_input => [
        target_intervals => {
            is => 'Text',
            gatk_param_name => '-targetIntervals',
            doc => 'The file of indels around which you wish to do realignment',
        },
        output_realigned_bam => {
            is_output => 1,
            is => 'Text',
            gatk_param_name => '-o',
            doc => 'The path to where you would like the realigned output bam',
        },
        input_bam => {
            is => 'Text',
            gatk_param_name => '-I',
            doc => 'The path to the original bam you would like to be realigned',
        },
        target_intervals_are_sorted => {
            is => 'Boolean',
            doc => 'If set to false, pass along --targetIntervalsAreNotSorted',
            default => 0,
        },
        reference_fasta => {
            is => 'Text',
            gatk_param_name => '-R',
            doc => "Reference Fasta" ,
        },
        index_bam => {
            is => 'Boolean',
            default => 1,
            doc => 'Index the bam after alignment.'
        },
        known => {
            is => 'Text',
            gatk_param_name => '--knownAlleles',
            doc => "Input VCF file(s) with known indels",
            is_optional => 1,
            is_many => 1,
        },
    ],
    has_optional_input => [
        _target_intervals_are_not_sorted => {
            is => 'Boolean',
            calculate_from => ['target_intervals_are_sorted'],
            calculate => q{ return !$target_intervals_are_sorted },
            gatk_param_name => '--targetIntervalsAreNotSorted',
        },
    ],
    has_param => [
        lsf_queue => {
            value => Genome::Config::get('lsf_queue_build_worker_alt'),
        },
    ],
};

sub help_brief {
    "Run GATK with the 'IndelRealigner' tool"
}

sub help_synopsis {
    return <<EOS
    gmt gatk indel-realigner --target-intervals some.bed --output-realigned-bam my_output_realigned.bam --input-bam my_existing.bam --reference-fasta my.fa
EOS
}

sub analysis_type {
    return 'IndelRealigner';
}

sub _shellcmd_extra_params {
    my $self = shift;

    unless ($self->target_intervals_are_sorted) {
        unless ($self->is_legacy_version($self->version)) {
            die $self->error_message("Version ".$self->version." does not support non-sorted target intervals.");
        }
    }

    return (
        input_files => [$self->known, $self->target_intervals, $self->input_bam, $self->reference_fasta],
        output_files => [$self->output_realigned_bam],
    );
}

sub _postprocess {
    my $self = shift;

    if ($self->index_bam) {
        my $aligned_bam = $self->output_realigned_bam;
        die $self->error_message("Couldn't find realigned bam at $aligned_bam!") unless -f $aligned_bam;

        my $rv = Genome::Model::Tools::Sam::IndexBam->execute(bam_file => $aligned_bam);
        die $self->error_message("Failed to run gmt sam index-bam on $aligned_bam") unless $rv->result == 1;
    }

    return 1;
}

sub _check_inputs {
    my $self = shift;

    unless ($self->target_intervals_are_sorted) {
        unless ($self->is_legacy_version($self->version)) {
            $self->error_message("Version ".$self->version." does not support non-sorted target intervals.");
            return;
        }
    }
    my @known = $self->known;
    if (@known) {
        for my $k (@known) {
            Genome::Sys->validate_file_for_reading(@known);
        }
    }
    Genome::Sys->validate_file_for_reading($self->target_intervals);
    Genome::Sys->validate_file_for_reading($self->input_bam);
    Genome::Sys->validate_file_for_reading($self->reference_fasta);
    Genome::Sys->validate_file_for_writing($self->output_realigned_bam);

    return 1;
}

1;
