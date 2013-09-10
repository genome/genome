#!/usr/bin/env genome-perl
use strict;
use warnings;
use Genome;

package Genome::Sys::Gateway::Command::Attach;

class Genome::Sys::Gateway::Command::Attach {
    is => 'Command::V2',
    has_input => [
        systems => { is => 'Genome::Sys::Gateway',
                    is_many => 1,
                    shell_args_position => 1,
                    doc => 'the system to attach'
                  },
    ],
    has_optional_param => [
        protocol  => { is => 'Text',
                      doc => 'override the protocol to be used to attach',
                    },
        rsync      => { is => 'Boolean',
                      default_value => 0,
                      doc => 'make a copy of the data for faster access'
                    },
    ],
    doc => 'attach the specified GMS gateways to the current GMS'
};

sub help_synopsis {
    return <<EOS
genome sys gateway attach GMS1
genome sys gateway attach GMS1 --protocol ftp
genome sys gateway attach GMS1 --protocol http
genome sys gateway attach GMS1 ABC123 XY1AB GMS2 --protocol nfs
genome sys gateway attach GMS1 ABC123 XY1AB GMS2 --protocol ssh
EOS
}

sub help_detail {
    return <<EOS
Attaches the specified other GMS gateway to the current one.

This means, if the system is named OTHER1, data will appear mounted at /opt/gms/OTHER1.

Using the "rsync" option will make a local copy of the attached data.  This can be 
quite large, but will improve performance during analysis.
EOS
}

sub execute {
    my $self = shift;
    my $protocol_override = $self->protocol;
    my @args = ($protocol_override ? ($protocol_override) : ());
    my $failures = 0;
    for my $sys ($self->systems) {
        $self->status_message("attaching remote GMS: " . $sys->id . "...");
        eval { $sys->attach(@args); };
        if ($@) {
            $self->error_message($@);
            $failures++;
        }
    }
    if ($self->rsync) {
        # this is terrible to not do in parallel if there are multiple systems
        for my $sys ($self->systems) {
            if (my $protocol = $sys->attached_via) {
                $sys->rsync($protocol);
            }
            else {
                $self->warning_message("skipping rsync of " . $self->__display_name__ . " because of errors attaching...\n");
            }
        }
    }
    return !$failures;
}

1;

