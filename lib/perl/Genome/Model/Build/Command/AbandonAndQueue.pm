package Genome::Model::Build::Command::AbandonAndQueue;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::AbandonAndQueue {
    is => 'Genome::Model::Build::Command::Abandon',
    has => [
        reason => {
            is => 'Text',
            doc => 'a brief note about why new builds are being enqueued',
        },
    ],
};

sub sub_command_sort_position { 6 }

sub help_brief {
    return "Abandon a build and queue the model for another build";
}

sub help_detail {
    return help_brief();
}

sub successfully_abandoned_callback {
    my $self = shift;
    my $build = shift;

    $build->model->build_requested(1, $self->reason);
    $self->status_message("Abandoned build (" . $build->__display_name__ . ") and queued model.");
}


1;
