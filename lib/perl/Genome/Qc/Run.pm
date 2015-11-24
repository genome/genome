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
        qc_genotype_vcf_file => {
            is => 'Genome::SoftwareResult::ImportedFile',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            default_value => Genome::Config::get('lsf_queue_build_worker_alt'),
        },
        lsf_resource => {
            default_value => Genome::Config::get('lsf_resource_qc_run'),
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
        qc_genotype_vcf_file => $self->qc_genotype_vcf_file,
    );
}

1;
