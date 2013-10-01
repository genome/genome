package Genome::File::Tsv;

use strict;
use warnings;
use Genome;

class Genome::File::Tsv {
    is => 'Genome::File::Xsv',
};

sub separator { "\t" }

1;

