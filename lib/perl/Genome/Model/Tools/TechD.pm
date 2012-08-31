package Genome::Model::Tools::TechD;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::TechD {
    is => ['Command'],
    is_abstract => 1,
};

sub help_detail {
    "These commands are used to perform Technology Development type anlaysis.";
}


1;
