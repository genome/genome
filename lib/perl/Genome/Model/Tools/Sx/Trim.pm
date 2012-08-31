package Genome::Model::Tools::Sx::Trim;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::Sx::Trim {
    is => 'Command',
    is_abstract => 1,
};

sub help_brief {
    return 'Trim sequences';
}

1;

