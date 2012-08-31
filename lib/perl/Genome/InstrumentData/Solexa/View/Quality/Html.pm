package Genome::InstrumentData::Solexa::View::Quality::Html;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Solexa::View::Quality::Html {
    is => 'Genome::View::Status::Html',
    has => [
        perspective => { value => 'quality' },
    ],
};

1;
