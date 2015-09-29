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
    my $reference_fasta_file = shift;

    unless ($self->aligner_class_name->prepare_annotation_index($self)) {
        $self->error_message("Failed to prepare annotation index.");
        return;
    }

    return $self;
}

sub _resolve_allocation_subdirectory_components {
    my $self = shift;

    return ('annotation_build_aligner_index_data',$self->reference_build->model->id,'reference_build'.$self->reference_build->id,'annotation_build'.$self->annotation_build_id);
}

1;
