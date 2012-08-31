package Genome::Model::Set::View::GoldSnpComparison::Html;

use strict;
use warnings;

use Genome;

class Genome::Model::Set::View::GoldSnpComparison::Html {
    is => 'UR::Object::View::Default::Html',
    has_constant => [
        perspective => {
            value => 'gold-snp-comparison',
        },
    ]
};

1;
