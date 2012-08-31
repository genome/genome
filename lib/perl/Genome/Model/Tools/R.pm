package Genome::Model::Tools::R;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::R {
    is => ['Command'],
};

sub help_brief {
    "A perl R-caller and associated R libraries."
}

sub help_detail {
    "This dir contains an R caller which uses the Statisics::R interface to open and close an R session. A temporary directory is created and passed to Statistics::R because this tool writes and deletes temporary files, which without a defined path, will end up in an operating system directory, and therefore could conflict if multiple calls are made simultaneously."
}

1;
