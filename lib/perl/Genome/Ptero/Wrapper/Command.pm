package Genome::Ptero::Wrapper::Command;

use strict;
use warnings FATAL => qw(all);
use Genome;

class Genome::Ptero::Wrapper::Command {
    is => 'Command::V2',
    has => [
        'command_class' => {
            is => 'Text',
            doc => 'class name of the Genome Command to execute',
        },
        'method' => {
            is => 'Text',
            valid_values => ['execute', 'shortcut'],
            doc => 'method to call on the Genome Command object',
        },
    ],

    doc => 'Ptero execution wrapper for Genome Command objects',
};


sub execute {
    die "This is not implemented yet, sorry";
}

1;
