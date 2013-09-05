#!/usr/bin/env genome-perl
use strict;
use warnings;
use Genome;

package Genome::Sys::Node::Command::Detach;

class Genome::Sys::Node::Command::Detach {
    is => 'Command::V2',
    has => [
        systems => { is => 'Genome::Sys::Node',
                is_many => 1,
                shell_args_position => 1,
                doc => 'the system to detach'
              },
    ],
    doc => 'detach the specified GMS gateways to the current GMS'
};

sub execute {
    my $self = shift;
    for my $sys ($self->systems) {
      $sys->detach();
    }
    return 1;
}

1;

