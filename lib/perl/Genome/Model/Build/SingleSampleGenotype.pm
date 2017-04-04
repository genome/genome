package Genome::Model::Build::SingleSampleGenotype;

use strict;
use warnings;

use Genome;

use Genome::Utility::Text;

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

sub _symlink_results {
    my $self = shift;
    my $destination = shift;

    my $top_directory = File::Spec->join($destination, $self->id);
    Genome::Sys->create_directory($top_directory);

    my $alignment = $self->merged_alignment_result;
    $self->_symlink_result($alignment, $top_directory, 'alignments');

    my $qc = $self->qc_result;
    $self->_symlink_result($qc, $top_directory, 'qc');

    $self->_symlink_haplotype_caller_results($top_directory);

    return $top_directory;
}

sub _symlink_result {
    my $self = shift;
    my $result = shift;
    my $location = shift;
    my $symlink_name = shift;

    my $symlink_location = File::Spec->join($location, $symlink_name);
    Genome::Sys->create_symlink($result->output_dir, $symlink_location);
}

sub _symlink_haplotype_caller_results {
    my $self = shift;
    my $top_directory = shift;

    my $variant_dir = File::Spec->join($top_directory, 'variants');
    Genome::Sys->create_directory($variant_dir);

    for my $hc_result ($self->haplotype_caller_result) {
        my $hc_basename = Genome::Utility::Text::sanitize_string_for_filesystem(
            join('-', $hc_result->intervals, $hc_result->id)
        );
        $self->_symlink_result($hc_result, $variant_dir, $hc_basename);
    }

    return $variant_dir;
}

sub _compare_output_files {
    my ($self, $other_build) = @_;

    my $temp = Genome::Sys->create_temp_directory();
    my $directory = $self->symlink_results($temp);
    my $other_directory = $other_build->symlink_results($temp);

    return $self->_compare_output_directories(
        $directory,
        $other_directory,
        $self->id,
        $other_build->id,
    );
}

sub _disk_usage_result_subclass_names {
    my $self = shift;

    my @disk_usage_result_classes = (
        'Genome::Qc::Result',
        'Genome::InstrumentData::AlignmentResult::Speedseq',
        'Genome::InstrumentData::AlignmentResult::Merged::Speedseq',
        'Genome::Model::SingleSampleGenotype::Result::HaplotypeCaller',
    );
    return \@disk_usage_result_classes;
}

1;
