package Genome::Model::Tools::Cmds;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Cmds {
    is => ['Command'],
};

sub help_brief {
    "Tools and R script to perform CMDS analysis."
}

sub help_detail {
    return <<EOS
Tools and R script to perform CMDS analysis.
EOS
}

1;
