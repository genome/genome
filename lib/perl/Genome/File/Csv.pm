package Genome::File::Csv;
use strict;
use warnings;
use Genome;

class Genome::File::Csv {
    is => 'Genome::File::Xsv',
    has_constant => [
        separator => { default_value => "," }
    ]
};

1;

