package Genome::Annotation::Expert::Fpkm::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::Annotation::Expert::Fpkm::Adaptor {
    is => "Genome::Annotation::AdaptorBase",

    has_planned_output => [
    ],
    has_output => [
        fpkm_file => {
            is => 'Path',
        },
    ],
};

sub name {
    'fpkm';
}

sub resolve_expert_specific_attributes_from_build {
    my $self = shift;
    $self->fpkm_file($self->resolve_fpkm_file);
    return;
}

sub resolve_fpkm_file {
    my $self = shift;
    return $self->fpkm_file_for_build($self->build_containing_fpkm);
}

sub build_containing_fpkm {
    my $self = shift;
    return $self->build->tumor_build;
}

sub fpkm_file_for_build {
    my ($self, $build) = @_;
    return File::Spec->join($build->data_directory, "expression", "genes.fpkm_tracking");
}

1;
