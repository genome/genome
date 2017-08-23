package Genome::InstrumentData::Command::AlignmentResult::Convert::Base;

use strict;
use warnings;

use Genome;
use Genome::Config;
use Genome::Sys::LSF::bsub qw();

class Genome::InstrumentData::Command::AlignmentResult::Convert::Base {
    is_abstract => 1,
    is => 'Command::V2',
    has_input => [
        alignment_result => {
            is => 'Genome::SoftwareResult',
            doc => 'result to convert',
        },
    ],
    doc => 'Convert an alignment result.',
};

sub _is_convertable_result {
    my $self = shift;

    my $alignment_result = $self->alignment_result;

    return unless $alignment_result->isa('Genome::InstrumentData::AlignedBamResult::Merged');
    return $alignment_result->can('filetype');
}

sub _lock {
    my $self = shift;

    my $result = $self->alignment_result;

    my $lock = Genome::Sys::LockProxy->new(
        resource => 'alignment-result-convert/' . $result->id,
        scope => 'site'
    )->lock(
        max_try => 1,
    );

    $self->fatal_message('Could not acquire lock for conversion') unless $lock;

    my ($input) = $result->inputs(name => 'filetype');
    if ($input) {
        #make sure we get any recent updates after acquiring the lock
        UR::Context->reload($input);
    }

    return $lock;
}

sub _is_currently_bam {
    my $self = shift;

    my $result = $self->alignment_result;

    my $filetype = $result->filetype;

    return 1 unless $filetype;
    return $filetype eq 'bam';
}

sub _is_currently_cram {
    my $self = shift;

    my $result = $self->alignment_result;

    my $filetype = $result->filetype;

    return 0 unless $filetype;
    return $filetype eq 'cram';
}

sub _run_conversion {
    my $self = shift;
    my $output_format_flag = shift;
    my $source_file = shift;
    my $destination_file = shift;

    my $reference = $self->result->reference_build;
    my $fasta = $reference->full_consensus_path('fa');

    my $cmd = ['view', '-T', $fasta, $output_format_flag, '-o', $destination_file, $source_file];

    my $guard = Genome::Config::set_env('lsb_sub_additional', 'docker(mgibio/samtools:1.3.1)');
    Genome::Sys::LSF::bsub::bsub(
        cmd => $cmd,
        queue => Genome::Config::get('lsf_queue_build_worker'),
        wait_for_completion => 1,
    );
}

sub _verify_file {
    my $self = shift;
    my $file = shift;

    return -e $file;
}

1;
