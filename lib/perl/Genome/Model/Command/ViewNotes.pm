package Genome::Model::Command::ViewNotes;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::ViewNotes {
    is => 'Genome::Notable::Command::ViewNotes',
    has => [
        notables => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            doc => 'notable objects on which to view the notes',
        },
    ],
    doc => 'view notes that have been set on models',
};

sub help_detail {
    return <<EOS
This command can be used to view the notes that have been added to a model.

For example this can be used to see why a build was requested for a model:
    genome model view-notes --note-type=build_requested <build_id>
EOS
}

1;
