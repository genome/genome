package Genome::DruggableGene::Command::Citation;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::Command::Citation {
    is => 'Command::Tree',
};

sub help_brief {
    "work with druggable gene citations"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 genome druggable-gene citation ...
EOS
}

sub help_detail {
    return <<EOS
A collection of commands to interact with druggable-gene citations.
EOS
}

1;
