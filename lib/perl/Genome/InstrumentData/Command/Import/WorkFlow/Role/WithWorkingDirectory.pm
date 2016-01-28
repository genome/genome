package Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory;

use strict;
use warnings;

use Genome;

role Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory {
    has => [
        working_directory => {
            is => 'Text',
            is_input => 1,
            doc => 'Working directory to put output files.',
        },
    ],
};

1;

