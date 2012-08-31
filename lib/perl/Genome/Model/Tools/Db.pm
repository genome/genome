package Genome::Model::Tools::Db;

use strict;
use warnings;

use Genome;     

class Genome::Model::Tools::Db {
    is => 'Command',
};

sub sub_command_sort_position { 14 }

sub help_brief {
    "database tools"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools db ...    
EOS
}

1;

