package Genome::Model::Tools::Analysis::454;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::454 {
    is => ['Genome::Model::Tools::Analysis'],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools for analysis of 454 data.",
}

1;

