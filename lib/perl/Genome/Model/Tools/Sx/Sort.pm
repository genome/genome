package Genome::Model::Tools::Sx::Sort;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::Sx::Sort {
    is  => 'Command::Tree',
    is_abstract => 1,
};

sub help_brief {
    return 'Sort sequences';
}

1;

