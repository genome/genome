package Genome::Model::Command::ViewNotes;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::ViewNotes {
    is => 'Command::V2',
    roles => [Genome::Role::Notable::Command::ViewNotes->create(notable_type => 'Genome::Model')],
    doc => 'view notes that have been set on models',
};

sub help_detail : Overrides(Genome::Role::Notable::Command::ViewNotes) {
    return <<EOS
This command can be used to view the notes that have been added to a model.

For example this can be used to see why a build was requested for a model:
    genome model view-notes --note-type=build_requested <build_id>
EOS
}

1;
