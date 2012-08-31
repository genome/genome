package Genome::Model::Build::Command::AddNote;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Command::AddNote {
    is => 'Genome::Notable::Command::AddNote',
    has => [
        notable => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'build to which a note will be added'
        },
    ],
};

sub help_detail {
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

