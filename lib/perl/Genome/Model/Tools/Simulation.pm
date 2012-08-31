package Genome::Model::Tools::Simulation;

use strict;
use warnings;

use Genome;            
our $VERSION = '0.01';

class Genome::Model::Tools::Simulation {
    is => 'Command',
    doc => 'A plethora of really useful stuff, you will want to check this one out ASAP.'
};

sub sub_command_sort_position { 16 }

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt simulation ...
EOS
}

1;

