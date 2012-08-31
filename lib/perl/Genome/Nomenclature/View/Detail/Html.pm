package Genome::Nomenclature::View::Detail::Html;

use strict;
use warnings;

use Genome;

class Genome::Nomenclature::View::Detail::Html {
    is => 'Genome::View::Status::Html',
    has => [
        perspective => { value => 'detail', },
    ],
};

1;
