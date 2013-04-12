package Genome::Config::CopyStrategy;

use warnings;
use strict;

use Genome;

class Genome::Config::CopyStrategy {
    is          => 'UR::Object',
    doc         => 'A copy strategy determines how config will be copied from the default location',
    is_abstract => 1,
};

sub copy_config {
    die('You must impletement copy_config!');
}

1;
