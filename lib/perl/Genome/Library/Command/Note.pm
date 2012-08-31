package Genome::Library::Command::Note;

use strict;
use warnings;

use Genome;

class Genome::Library::Command::Note {
    #is => 'UR::Object::Command::Note',
    is => 'Command',
    doc => 'notes for libraries',
};

sub sub_command_sort_position { 5 }

1;

