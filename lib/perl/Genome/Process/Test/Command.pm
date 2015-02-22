package Genome::Process::Test::Command;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Process::Test::Command {
    is => ['Command::V2'],
    has => {
        statement => {
            is => 'Text',
            is_input => 1,
        },
        some_output => {
            is => 'Text',
            is_output => 1,
        }
    },
};

sub execute {
    my $self = shift;
    $self->error_message("STATEMENT: %s", $self->statement);
    $self->some_output($self->statement);
    return 1;
}

1;
