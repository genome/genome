package Genome::Model::ImportedVariationList::Command::BuildForReference;

use strict;
use warnings;
use Genome;

class Genome::Model::ImportedVariationList::Command::BuildForReference {
    is => 'Command::V2',
    has => [
        reference_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            shell_args_position => 1,
            doc => 'Reference build for which we\'re finding a variation build',
        },
    ],
    doc => 'Finds a variation build that goes with the supplied reference build',
};

sub help_detail {
    return 'Finds a variation build that goes with the supplied reference build';
}

sub execute {
    my $self = shift;

    my $build = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($self->reference_build);

    if ($build) {
        $self->status_message("Found imported variation build " . $build->__display_name__ . 
            " which should be used in conjunction with supplied reference " . $self->reference_build->__display_name__);
    }
    else {
        $self->status_message("Found no imported variation build for reference " . $self->reference_build->__display_name__);
    }

    return 1;
}

1;

