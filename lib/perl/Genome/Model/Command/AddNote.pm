package Genome::Model::Command::AddNote;

use strict;
use warnings;
use Genome;

class Genome::Model::Command::AddNote {
    is => 'Command::V2',
    roles => [Genome::Role::Notable::Command::AddNote->create(notable_type => 'Genome::Model')],
};

sub help_detail : Overrides(Genome::Role::Notable::Command::AddNote) {
    return <<EOS
add a note to the specified model
EOS
}

sub help_brief {
    return help_detail();
}

sub help_synopsis {
    return <<EOS
genome model add-note 123456 --header-text 'header' --body-text 'body'
EOS
}

1;

