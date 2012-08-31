package Genome::Model::Tools::Relationship;

use strict;
use warnings;

use Genome;            
our $VERSION = '0.01';

class Genome::Model::Tools::Relationship {
    is => 'Command',
    doc => 'Analysis you might want to perform on related individuals'
};

sub sub_command_sort_position { 16 }

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt relationship ...
EOS
}

1;

