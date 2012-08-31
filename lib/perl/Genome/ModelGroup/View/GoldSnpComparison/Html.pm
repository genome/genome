package Genome::ModelGroup::View::GoldSnpComparison::Html;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::View::GoldSnpComparison::Html {
    is => 'UR::Object::View::Default::Html',
    has_constant => [
        perspective => {
            value => 'gold-snp-comparison',
        },
    ]
};

1;
