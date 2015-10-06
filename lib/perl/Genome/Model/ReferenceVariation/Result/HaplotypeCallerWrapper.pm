package Genome::Model::ReferenceVariation::Result::HaplotypeCallerWrapper;

class Genome::Model::ReferenceVariation::Result::HaplotypeCallerWrapper {
    is => 'Genome::Command::DelegatesToResult',
    has_input => {
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult',
            doc => 'The alignments on which to run the haplotype caller',
        },
        intervals => {
            is_many => 1,
            is => 'Text',
            doc => 'regions to restrict this result',
        },
        gvcf_gq_bands => {
            is_many => 1,
            is => 'Number',
            doc => 'GQ thresholds for this result',
        },
        emit_reference_confidence => {
            is => 'Text',
            doc => 'ERC option for this result',
        },
        haplotype_caller_version => {
            is => 'Text',
            doc => 'Version of GATK to use',
        },
    },
};

sub result_class { return 'Genome::Model::ReferenceVariation::Result::HaplotypeCaller'; }

1;
