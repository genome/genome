package Genome::Qc::Run;

use strict;
use warnings;
use Genome;

class Genome::Qc::Run {
    is => 'Genome::Command::DelegatesToResult',
    has_input => [
        config_name => {
            is => 'Text',
        },
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
    ],
    has_param => [
        lsf_queue => {
            default_value => Genome::Config::get('lsf_queue_build_worker_alt'),
        },
        lsf_resource => {
            value => &bsub_rusage,
        },
    ],
};

sub bsub_rusage {
    return sprintf("-q %s -R 'span[hosts=1] select[mem>16000] rusage[mem=16000]' -M 16000000 -n 4", Genome::Config::get('lsf_queue_build_worker_alt')),
}

sub result_class {
    return "Genome::Qc::Result";
}

sub input_hash {
    my $self = shift;
    return (
        alignment_result => $self->alignment_result,
        config_name => $self->config_name,
    );
}

1;
