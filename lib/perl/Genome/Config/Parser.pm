package Genome::Config::Parser;

use warnings;
use strict;

class Genome::Config::Parser {
    is          => 'UR::Object',
    is_abstract => 1,
    doc         => 'interface for config file parsers to follow',
};

sub parse { die('You must implement a parse method!'); }

1;
