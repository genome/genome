package Genome::Model::Build::ReferenceAlignment::View::Coverage::Html;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::ReferenceAlignment::View::Coverage::Html {
    is => 'Genome::View::Status::Html',
    has => [
        perspective => { value => 'coverage' },
    ],
};

1;
