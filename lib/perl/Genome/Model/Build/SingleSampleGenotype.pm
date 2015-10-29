package Genome::Model::Build::SingleSampleGenotype;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::SingleSampleGenotype {
    is => ['Genome::Model::Build'],
    has_output => [
        merged_alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult::Merged',
            doc => 'the result of running alignment',
        },
        qc_result => {
            is => 'Genome::Qc::Result',
            doc => 'the result of qc analysis on the alignment',
        },
        haplotype_caller_result => {
           is => 'Genome::Model::SingleSampleGenotype::Result::HaplotypeCaller',
           doc => 'the result of running the haplotype caller',
           is_many => 1,
        },
        buckets_result => {
            is => 'Genome::Model::Build::ReferenceSequence::Buckets',
            doc => 'the buckets used for parallelization',
        },
    ],
};

sub workflow_name {
    my $self = shift;
    return $self->build_id . ' Single Sample Genotype Pipeline';
}

sub calculate_estimated_kb_usage {
    my $self = shift;

    return 50_000; #just need space for logs
}

sub reference_being_replaced_for_input {
    my ($self, $input) = @_;

    if($input->name eq "target_region_set"){
        return 1; #not used
    }

    if($input->name eq "region_of_interest_set"){
        my $rsb = $self->reference_sequence_build;
        my $roi_reference = $input->value->reference;
        unless ($roi_reference) {
            return;
        }

        if ($roi_reference and !$rsb->is_compatible_with($roi_reference)) {
            my $converter_exists =  Genome::Model::Build::ReferenceSequence::Converter->exists_for_references(
                $roi_reference, $rsb,
            );

            if ($converter_exists) {
                return 1;
            }
        }
    }

    if($input->name eq 'instrument_data') {
        return 1; #we'll re-align
    }

    return;
}

1;
