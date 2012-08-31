package Genome::Model::Tools::Dbsnp;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dbsnp {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for importing dbSnp files.'
} 

sub help_synopsis {
    return <<"EOS"
gmt dbsnp ...
EOS
}

1;
