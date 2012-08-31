package Genome::Model::Tools::Sx::Bin;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Bin {
    is => 'Command::Tree',
    is_abstract => 1,
};

sub help_brief {
    return 'Bin sequences';
}

1;

