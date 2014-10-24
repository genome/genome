package Genome::Model::Command::AddNote;

use strict;
use warnings;
use Genome;

class Genome::Model::Command::AddNote {
    is => 'Genome::Notable::Command::AddNote',
    has => [
        notable => {
            is => 'Genome::Model',
            shell_args_position => 1,
            doc => 'build to which a note will be added'
        },
    ],
};

sub help_detail {
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

