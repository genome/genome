package Genome::Model::Build::Command::ViewNotes;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::ViewNotes {
    is => 'Genome::Notable::Command::ViewNotes',
    has => [
        notables => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'notable objects on which to view the notes',
        },
    ],
    doc => 'view notes that have been set on builds',
};

sub help_detail {
    return <<EOS
This command can be used to view the notes that have been added to a build.

For example this can be used to see why a build was not startable:
    genome model build view-notes --note-type=Unstartable <build_id>
EOS
}

1;