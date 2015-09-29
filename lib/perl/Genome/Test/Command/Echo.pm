package Genome::Test::Command::Echo;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Test::Command::Echo {
    is => 'Command::V2',

    has_input => [
        input_statement => {
            is => 'Text',
        },
        fail_intentionally => {
            is => 'Boolean',
            default => 0,
        },
    ],
    has_optional_output => [
        output_statement => {
            is => 'Text',
        },
    ],
    doc => 'Echos a statement',
};

sub shortcut {
    my $self = shift;
    return $self->execute();
}

sub execute {
    my $self = shift;

    die "I was told to fail, so here goes!" if $self->fail_intentionally;

    $self->status_message("I'm about to echo a statement");
    $self->status_message($self->input_statement);
    $self->status_message("I just echoed a statement, did you like it?");

    $self->output_statement($self->input_statement);

    return 1;
}

1;
