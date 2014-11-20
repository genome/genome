package Genome::Process::Test::Process;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Process::Test::Process {
    is => ['Genome::Process'],
    has => [
        statement => {
            is => 'Text',
            is_input => 1,
        },
    ],
};

1;
