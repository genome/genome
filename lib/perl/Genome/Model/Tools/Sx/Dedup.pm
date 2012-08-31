package Genome::Model::Tools::Sx::Dedup;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Dedup {
    is  => 'Command::Tree',
};

sub help_brief {
    return 'Deduplicate sequences';
}

1;

