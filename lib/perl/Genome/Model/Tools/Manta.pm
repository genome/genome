package Genome::Model::Tools::Manta;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Manta {
    is  => 'Command::Tree',
    doc => 'Toolkit for executing the Manta structural variant call',
};

sub help_detail {
    "These commands are used to execute Manta and related processes";
}
