package Genome::Model::Tools::Sx::Split;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Split {
    is  => 'Command::Tree',
    is_abstract => 1,
};

sub help_brief {
    return 'split sequences';
}

1;

