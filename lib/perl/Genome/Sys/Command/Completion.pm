package Genome::Sys::Command::Completion;

use Genome;
use strict;
use warnings;

class Genome::Sys::Command::Completion {
    is => 'Genome::Command::Base',
    is_abstract => 1,
    doc => 'Commands to assist with shell tab completion.',
};

1;
