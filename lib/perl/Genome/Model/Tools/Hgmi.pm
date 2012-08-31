package Genome::Model::Tools::Hgmi;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Hgmi {
    is => ['Command'],
};

sub help_brief {
    "tools to work for the Hgmi annotation pipeline"
}


sub sub_command_sort_position { 16 }

sub help_detail {
    return <<EOS
need to fill out the hgmi help detail
EOS
}


1;

