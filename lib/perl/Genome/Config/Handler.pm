package Genome::Config::Handler;

use warnings;
use strict;

class Genome::Config::Handler {
    is          => 'UR::Object',
    doc         => 'interface for config handlers to follow',
    is_abstract => 1,
};

sub valid_params { die('must implement valid_params!'); }

sub get_config { die('must implement get_config!'); }

1;
