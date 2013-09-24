package Genome::Model::Tools::ChimeraScan;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ChimeraScan {
    is  => 'Command::Tree',
    doc => 'Toolkit for ChimeraScan related process',
};

sub help_detail {
    "These commands are used to perform ChimeraScan related process";
}

1;
