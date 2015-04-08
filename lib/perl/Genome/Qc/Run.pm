package Genome::Qc::Run;

use strict;
use warnings;
use Genome;

class Genome::Qc::Run {
    is => 'Genome::Command::DelegatesToResult',
    has => [
        config_name => {
            is => 'Text',
        },
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
    ],
    has_param => [
        lsf_resource => {
            value => q{-R 'select[mem>16000] rusage[mem=16000]' -M 16000000},
        },
    ],
};

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
