package Genome::File::Csv;
use strict;
use warnings;
use Genome;

class Genome::File::Csv {
    is => 'Genome::File::Xsv',
};

sub separator { ',' }

1;

