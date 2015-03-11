package Genome::Qc::Run;

use strict;
use warnings;
use Genome;

class Genome::Qc::Run {
    is => 'Genome::Command::DelegatesToResult',
    has => [
        #qc_config => {
        #    is => 'Genome::Qc::Config',
        #},
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
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
        #qc_config => $self->qc_config, #enable when db-backed qc_config is available
    );
}
