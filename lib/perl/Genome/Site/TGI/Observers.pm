package Genome::Site::TGI::Observers;

use strict;
use warnings;

# NOTE: the global syslog observers are in Genome::Sys::Log module.

UR::Observer->register_callback(
    subject_class_name => 'UR::Object::Type',
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

my $unarchive_observer;
$unarchive_observer = UR::Observer->register_callback(
    subject_class_name => 'UR::Object::Type',
    aspect => 'load',
    callback => sub {
        my $meta = shift;
        my $class_name = $meta->class_name;
        if ($class_name eq 'Genome::InstrumentData::Command::Unarchive') {
            require Genome::Site::TGI::Extension::UnarchiveInstrumentData;
            UR::Observer->unregister_callback(id => $unarchive_observer);
        }
    },
);

my $cqid_observer;
$cqid_observer = UR::Observer->register_callback(
    subject_class_name => 'UR::Object::Type',
    aspect => 'load',
    callback => sub {
        my $meta = shift;
        my $class_name = $meta->class_name;
        if ($class_name eq 'Genome::Config::Command::ConfigureQueuedInstrumentData') {
            require Genome::Site::TGI::Extension::ConfigureQueuedInstrumentData;
            UR::Observer->unregister_callback(id => $cqid_observer);
        }
    },
);

1;

