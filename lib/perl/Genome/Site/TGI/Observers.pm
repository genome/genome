package Genome::Site::TGI::Observers;

use strict;
use warnings;

# NOTE: the global syslog observers are in Genome::Sys::Log module.

UR::Object::Type->add_observer(
    aspect => 'load',
    callback => sub {
        my $meta = shift;
        my $class_name = $meta->class_name;
        if ($class_name eq 'Genome::ModelGroup') {
            require Genome::Site::TGI::Observers::ModelGroup;
        } elsif ($class_name eq 'Genome::Project') {
            require Genome::Site::TGI::Observers::Project;
        }
        die $@ if $@;
    },
);

1;

