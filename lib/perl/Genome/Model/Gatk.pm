package Genome::Model::Tools::Gatk;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gatk {
    is => 'Command::Tree',
};

sub help_detail {
    return 'Commmands fo using the GATK';
}

1;
