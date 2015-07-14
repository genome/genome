package Genome::Model::Build::ReferenceSequence::AnnotationIndex;

use Genome;
use warnings;
use strict;


class Genome::Model::Build::ReferenceSequence::AnnotationIndex {
    is => ['Genome::Model::Build::ReferenceSequence::IndexBase'],
    has => [
        annotation_build         => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'annotation_build_id',
        },
        annotation_name         => { via => 'annotation_build', to => 'name', is_mutable => 0, is_optional => 1 },
    ],
    has_input => [
        annotation_build_id      => {
            is => 'Number',
            doc => 'the annotation to use by id',
        },
    ],
};

sub _working_dir_prefix {
    "annotation-index";
}

sub __display_name__ {
    my $self = shift;
    my @class_name = split("::", $self->class);
    my $class_name = $class_name[-1];

    return sprintf("%s for reference build %s and annotation build %s with %s, version %s, params='%s'",
                   $class_name,
                   $self->reference_name,
                   $self->annotation_name,
                   $self->aligner_name,
                   $self->aligner_version,
                   $self->aligner_params || ""
               );
}

sub create {
    my $class = shift;
    my %p = @_;

    my $aligner_class = 'Genome::InstrumentData::AlignmentResult::'  . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($p{aligner_name});

    $class->debug_message(sprintf("Resolved aligner class %s, making sure it's real and can be loaded.", $aligner_class));
    unless ($aligner_class->class) {
        $class->error_message(sprintf("Failed to load aligner class (%s).", $aligner_class));
        return;
    }

    my $self = $class->SUPER::create(%p);
    return unless $self;
    $self->aligner_class_name($aligner_class);

    $self->debug_message("Prepare staging directories...");
    unless ($self->_prepare_staging_directory) {
        $self->error_message("Failed to prepare working directory");
        return;
    }

    unless ($self->_prepare_index) {
        $self->error_message("Failed to prepare annotation index!");
        return;
    }

    unless ($self->generate_dependencies_as_needed()) {
        $self->error_message("Failed to create AnnotationIndex objects for dependencies");
        return;
    }

    return $self;
}

sub generate_dependencies_as_needed {
    my $self = shift;

    # if the reference is a compound reference
    if ($self->reference_build->append_to and $self->_supports_multiple_reference) {
        die('Compound references are not currently supported in '. __PACKAGE__);
    }

    return 1;
}

sub _prepare_index {
    my $self = shift;

    my $reference_fasta_file;
    if ($self->_supports_multiple_reference) {
        $reference_fasta_file = $self->reference_build->primary_consensus_path('fa');
    } else {
        $reference_fasta_file = $self->reference_build->full_consensus_path('fa');
    }

    unless (-s $reference_fasta_file) {
        $self->error_message(sprintf("Reference fasta file %s does not exist", $reference_fasta_file));
        return;
    }

    $self->debug_message(sprintf("Confirmed non-zero reference fasta file is %s", $reference_fasta_file));

    unless ($self->aligner_class_name->prepare_annotation_index($self)) {
        $self->error_message("Failed to prepare annotation index.");
        return;
    }

    my $output_dir = $self->output_dir || $self->_prepare_output_directory;
    $self->debug_message("Alignment output path is $output_dir");

    unless ($self->_promote_data)  {
        $self->error_message("Failed to de-stage data into output path " . $self->output_dir);
        return;
    }

    $self->_reallocate_disk_allocation;

    $self->debug_message("Prepared alignment reference index!");

    return $self;
}

sub _resolve_allocation_subdirectory_components {
    my $self = shift;

    return ('annotation_build_aligner_index_data',$self->reference_build->model->id,'reference_build'.$self->reference_build->id,'annotation_build'.$self->annotation_build_id);
}

1;
