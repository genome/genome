package Genome::Model::Build::Command::AddNote;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Command::AddNote {
    is => 'Command::V2',
    roles => [Genome::Role::Notable::Command::AddNote->create(notable_type => 'Genome::Model::Build')],
};

sub help_detail : Overrides(Genome::Role::Notable::Command::AddNote) {
    return <<EOS
This command adds a note to the specified build
EOS
}

sub help_brief {
    return help_detail();
}

sub help_synopsis {
    return <<EOS
genome model build add-note 123456 --header-text 'header' --body-text 'body'
EOS
}

1;

