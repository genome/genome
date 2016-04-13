package Genome::Model::Build::SingleSampleGenotype;

use strict;
use warnings;

use Genome;
use Set::Scalar;
use Data::Compare;

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

sub _compare_output_files {
    my ($self, $other_build) = @_;

    my %diffs = (
        $self->_compare_merged_alignment_result($other_build),
        $self->_compare_haplotype_caller_results($other_build),
        $self->_compare_qc_metrics($other_build),
    );

    return %diffs;
}

sub _compare_merged_alignment_result {
    my ($self, $other_build) = @_;

    my $blessed_merged_alignment_result = $self->merged_alignment_result->bam_file;
    my $other_merged_alignment_result = $other_build->merged_alignment_result->bam_file;
    my $diff_result = !$self->file_comparer->compare($blessed_merged_alignment_result, $other_merged_alignment_result);

    unless ($diff_result) {
        return ('merged_alignment_result' => "files are not the same (diff -u $blessed_merged_alignment_result $other_merged_alignment_result)");
    }

    return;
}

sub _compare_haplotype_caller_results {
    my ($self, $other_build) = @_;

    my @blessed_haplotype_caller_results = $self->haplotype_caller_result;
    my @other_haplotype_caller_results = $other_build->haplotype_caller_result;

    my %diffs;
    for my $blessed_haplotype_caller_result (@blessed_haplotype_caller_results) {
        my $blessed_haplotype_caller_result_intervals = Set::Scalar->new($blessed_haplotype_caller_result->intervals);
        my ($other_haplotype_caller_result) = grep {
            my $other_haplotype_caller_result_intervals = Set::Scalar->new($_->intervals);
            $other_haplotype_caller_result_intervals->is_equal($blessed_haplotype_caller_result_intervals);
        } @other_haplotype_caller_results;
        my $diff_result = !$self->file_comparer->compare($blessed_haplotype_caller_result, $other_haplotype_caller_result);

        unless ($diff_result) {
            $diffs{sprintf("haplotype_caller_result-%s", $blessed_haplotype_caller_result->id)} =
                "files are not the same (diff -u $blessed_haplotype_caller_result->vcf_file $other_haplotype_caller_result->vcf_file)";
        }
    }

    return %diffs;
}

sub _compare_qc_metrics {
    my ($self, $other_build) = @_;

    my %blessed_qc_metrics = $self->qc_result->get_metrics;
    my %other_qc_metrics = $other_build->qc_result->get_metrics;

    unless (Compare(\%blessed_qc_metrics, \%other_qc_metrics)) {
        return ('qc_result' => "qc metrics are not the same");
    }

    return;
}

1;
