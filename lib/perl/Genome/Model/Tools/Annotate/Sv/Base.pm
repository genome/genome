package Genome::Model::Tools::Annotate::Sv::Base;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Annotate::Sv::Base{
    is => "UR::Object",
    has => [
    ],
};

sub process_breakpoint_list {
    #override in subclass;
    #interface for sv annotators
};

1;
