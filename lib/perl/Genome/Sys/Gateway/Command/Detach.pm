#!/usr/bin/env genome-perl

package Genome::Sys::Gateway::Command::Detach;
use strict;
use warnings;
use Genome;

class Genome::Sys::Gateway::Command::Detach {
    is => 'Command::V2',
    has_input => [
        systems => { is => 'Genome::Sys::Gateway',
                    is_many => 1,
                    shell_args_position => 1,
                    doc => 'the system to detach'
                  },
    ],
    has_optional_param => [
        protocol  => { is => 'Text',
                      doc => 'override the protocol to be used to detach',
                    },
    ],
    doc => 'detach the specified GMS gateways to the current GMS'
};

sub help_synopsis {
    return <<EOS
genome sys node detach GMS1
genome sys node detach GMS1 --protocol ftp
genome sys node detach GMS1 --protocol http
genome sys node detach GMS1 ABC123 XY1AB GMS2 --protocol nfs
genome sys node detach GMS1 ABC123 XY1AB GMS2 --protocol ssh
EOS
}

sub help_detail {
    return <<EOS
Detaches the specified other GMS gateway from the current one.

This means, if the system is named OTHER1, data mounted at /opt/gms/OTHER1 will become unavailable.
EOS
}

sub execute {
    my $self = shift;
    my $protocol_override = $self->protocol;
    my @args = ($protocol_override ? ($protocol_override) : ());
    my $failures = 0;
    for my $sys ($self->systems) {
        $self->status_message("detaching remote GMS: " . $sys->id . "...");
        eval { $sys->detach(@args); };
        if ($@) {
            $self->error_message($@);
            $failures++;
        }
    }
    return !$failures;
}

1;

