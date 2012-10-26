package Genome::DataSource::LocalDataSource;

use strict;
use warnings;

use Carp qw(croak);
use File::Temp;

class Genome::DataSource::LocalDataSource {
    is => [ 'UR::DataSource::SQLite', 'UR::Singleton' ],
};

sub known_classes {
    return qw(
        Genome::Disk::Allocation
        Genome::Disk::Assignment
        Genome::Disk::Group
        Genome::Disk::Volume
    );
}

our @hijacked_classes;
sub import {
    my $class = shift;
    my @classes_to_hijack = @_;

    my $self = $class->get();
    for my $class_to_hijack (@classes_to_hijack) {
        $self->hijack_class($class_to_hijack) or die;
    }
}

sub server {
    our $PID;

    # Store file in environment variable so that it can persist to subshells
    # which is needed since allocation create subshells to perform _create,
    # etc.
    unless ($ENV{DS_LDS_FILE}) {
        $PID = $$;
        (undef, $ENV{DS_LDS_FILE}) = File::Temp::tempfile(
            'genome_disk_localdatasource_XXXXXX',
            DIR    => $ENV{GENOME_TEST_TEMP},
            OPEN   => 0,
            UNKINK => 0, # UNLINK Fix in UR too!?
            TMPDIR => 1,
            SUFFIX => '.sqlite3',
        );
        unless ($ENV{DS_LDS_FILE}) {
            Carp::croak(sprintf('Failed to create tempfile for %s!', __PACKAGE__));
        }
    }
    return $ENV{DS_LDS_FILE};
}

sub _dbi_connect_args {
    my $self = shift;

    # Have to override _dbi_connect_args to enable adequate locking for concurrency.

    my @connection;
    $connection[0] = $self->dbi_data_source_name;
    $connection[1] = $self->login;
    $connection[2] = $self->auth;
    $connection[3] = { RaiseError => 0, sqlite_use_immediate_transaction => 1 };

    return @connection;
}

# Don't print warnings about loading up the DB if running in the test harness
sub _dont_emit_initializing_messages {
    my($dsobj, $msg) = @_;

    if ($msg =~ m/^Re-creating/) {
        $_[1] = undef;
    }
}

if ($ENV{'HARNESS_ACTIVE'}) {
    # don't emit messages while running in the test harness
    __PACKAGE__->warning_messages_callback(\&_dont_emit_initializing_messages);
}

sub class_to_filename {
    my $class = shift;
    my $filename = $class;
    $filename =~ s/::/\//g;
    $filename .= '.pm';
    return $filename;
}

sub hijack_class {
    my $self  = shift;
    my $class = shift;

    my $filename = class_to_filename($class);
    if ($INC{$filename}) {
        croak "$class has already been loaded!";
    }

    $class->class();
    unless ($INC{$filename}) {
        die "$class did not load from expected filename ($filename) hijack cannot be reliable!";
    }

    my $meta = $class->__meta__;
    my $old_ds_id = $meta->data_source_id;
    $meta->data_source_id($self->id);
    if ($meta->data_source_id eq $old_ds_id) {
        die "$class data_source not changed!";
    }

    return 1;
}

END {
    our $PID;
    if ($PID && $PID == $$ && $ENV{DS_LDS_FILE} && -e $ENV{DS_LDS_FILE}) {
        unlink $ENV{DS_LDS_FILE};
    }
}


1;
