package Genome::Show;

use strict;
use warnings;


class Genome::Show {
    is => 'Command::V2',

    has => [
        target => {
            is => 'Text',
            shell_args_position => 1,
            doc => 'The string, id, or partial id to search for.',
        },
        color => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display report in color.'
        },
    ],
};


1;
