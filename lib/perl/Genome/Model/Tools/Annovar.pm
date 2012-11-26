package Genome::Model::Tools::Annovar;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Annovar {
    is => 'Command::Tree',
};

sub help_detail {
    return 'Annotate with ANNOVAR';
}

1;
