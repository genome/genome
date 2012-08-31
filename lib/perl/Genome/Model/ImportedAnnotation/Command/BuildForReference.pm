package Genome::Model::ImportedAnnotation::Command::BuildForReference;

use strict;
use warnings;
use Genome;

class Genome::Model::ImportedAnnotation::Command::BuildForReference {
    is => 'Command::V2',
    has => [
        reference_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            shell_args_position => 1,
            doc => 'Reference sequence build for which we\'re finding a compatible annotation build',
        },
    ],
    doc => 'Finds a compatible annotation reference build for a provided reference build',
};

sub help_detail {
    return 'Finds a compatible annotation reference build for a provided reference build';
}

sub execute {
    my $self = shift;

    my $build = Genome::Model::ImportedAnnotation->annotation_build_for_reference($self->reference_build);

    if ($build) {
        $self->status_message("Found imported annotation build " . $build->__display_name__ . 
            " which should be used in conjunction with supplied reference " . $self->reference_build->__display_name__);
    }
    else {
        $self->status_message("Found no imported annotaiton build for reference " . $self->reference_build->__display_name__);
    }

    return 1;
}

1;

