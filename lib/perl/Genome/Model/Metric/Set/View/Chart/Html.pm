package Genome::Model::Metric::Set::View::Chart::Html;

use strict;
use warnings;

use Genome;

class Genome::Model::Metric::Set::View::Chart::Html {
    is => 'Genome::View::Status::Html',
    has => [
        perspective => { value => 'chart', },
    ],
};

1;
