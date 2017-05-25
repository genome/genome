package Genome::InstrumentData::Command::AlignmentResult::Import;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::AlignmentResult::Import {
    is => 'Command::V2',
    has_input => [
        alignment_file => {
            is => 'FilePath',
            doc => 'The BAM/CRAM file of alignments to import'
        },
        sample => {
            is => 'Genome::Sample',
            doc => 'The sample from which the data for these alignments were taken',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference sequence build to which the alignments were made',
        },
    ],
    has_optional_input => [
        description => {
            is => 'Text',
            doc => 'brief information about the alignment data',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'Specific instrument data used for these alignments (if any)',
            is_many => 1,
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'annotation build used for alignments (if any)',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->_verify_alignment_file;

    my $result_params = $self->_collect_result_params;
    my $result = Genome::InstrumentData::AlignmentResult::Merged::External->create(
        %$result_params
    );

    return $result;
}

sub _verify_alignment_file {
    my $self = shift;

    my $alignment_file = $self->alignment_file;
    unless (-e $alignment_file) {
        $self->fatal_message('Error reading %s: No such file or directory.', $alignment_file);
    }

    return 1;
}

sub _collect_result_params {
    my $self = shift;

    my %result_params = (
            original_file_path => $self->alignment_file,
            sample => $self->sample,
            reference_build => $self->reference_build,
    );

    if (my $description = $self->description) {
        $result_params{description} = $description;
    }

    if (my @instrument_data = $self->instrument_data) {
        $result_params{instrument_data} = \@instrument_data;
    }

    if (my $annotation_build = $self->annotation_build) {
        $result_params{annotation_build} = $annotation_build;
    }

    return \%result_params;
}

1;
