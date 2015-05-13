#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 1;

use UR;

my $created_observers = 0;
UR::Observer->add_observer(
    subject_class_name => 'UR::Observer',
    aspect => 'create',
    callback => sub {
        my ($observer, $aspect) = @_;

        $created_observers++;

        my $object;
        if (not $observer->subject_class_name) {
            $object = 'all classes';
        }
        elsif (not $observer->subject_id) {
            $object = 'class ' . $observer->subject_class_name;
        }
        else {
            $object = sprintf('%s object (ID: %s)', $observer->subject_class_name, $observer->subject_id);
        }

        my $message;
        if (not $observer->aspect) {
            $message = sprintf(q(Observer created on %s.), $object);
        }
        else {
            $message = sprintf(q(Observer created on %s watching '%s' aspect.), $object, $observer->aspect);
        }

        diag $message;
    },
);

require Genome::Model::Command::Services::WebApp::Core;
Genome::Model::Command::Services::WebApp::Core->import();

# If an observer is created during import of Genome packages it should use the
# `UR::Observer->register_callback` API so that they are rollback-safe.
is($created_observers, 0, 'no observers should be created during import of Genome packages');
