package Genome::File::Tsv;

use strict;
use warnings;
use Genome;

class Genome::File::Tsv {
    is => 'Genome::File::Base',
    has_constant => [
        separator => { default_value => "\t" }
    ]
};

1;

