package Genome::Model::Tools::Dindel;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dindel {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for running the hacky convergence of python and c known as Dindel. Bet you thought I was gonna say Tophat.'
} 

sub help_synopsis {
    return <<"EOS"
gmt dbsnp ...
EOS
}

1;
