package Genome::Annotation::Adaptor::SomaticVariation;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Adaptor::SomaticVariation {
    is => 'Genome::Annotation::Adaptor',
    has_input => [
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
        },
    ],
};

sub resolve_bam_results {
    my $self = shift;
    my @bam_results;
    for my $type qw(normal_build tumor_build) {
        push @bam_results, $self->build->$type->merged_alignment_result;
    }
    return \@bam_results;
}

sub resolve_annotation_build {
    my $self = shift;
    return $self->build->annotation_build;
}

1;

