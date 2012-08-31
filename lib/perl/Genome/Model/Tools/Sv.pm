package Genome::Model::Tools::Sv;

use strict;
use warnings;

use Genome;            
our $VERSION = '0.01';

class Genome::Model::Tools::Sv {
    is => 'Command',
    doc => 'Tools for working with SV files of various kinds.'
};

sub sub_command_sort_position { 16 }

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt sv ...
EOS
}

1;

