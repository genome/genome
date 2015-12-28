package Genome::Model::Build::Command::ViewNotes;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::ViewNotes {
    is => 'Command::V2',
    roles => [Genome::Role::Notable::Command::ViewNotes->create(notable_type => 'Genome::Model::Build')],
    doc => 'view notes that have been set on builds',
};

sub help_detail : Overrides(Genome::Role::Notable::Command::ViewNotes) {
    return <<EOS
This command can be used to view the notes that have been added to a build.

For example this can be used to see why a build was not startable:
    genome model build view-notes --note-type=Unstartable <build_id>
EOS
}

1;
