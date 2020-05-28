package Genome::Disk::Command::Volume::List;

use strict;
use warnings;

use Genome;
use IO::File;

class Genome::Disk::Command::Volume::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Disk::Volume',
        },
        show => { 
            default_value => 'mount_path,disk_group_names,total_kb,percent_used,percent_allocated', 
        },
        accurate_size => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Mount each displayed volume to display used/free sizes accurately',
        },
    ],
    doc => 'Lists Genome::Disk::Volume objects',
};

sub execute {
    my $self = shift;

    require Genome::Disk::Group;
    require Genome::Disk::Volume;
    no warnings 'redefine';
    local *Genome::Disk::Group::validate = sub { 1 };  # For just listing, don't validate them
    local *Genome::Disk::Volume::used_kb = $self->make_is_mounted_wrapper(\&Genome::Disk::Volume::used_kb);
    local *Genome::Disk::Volume::percent_used = $self->make_is_mounted_wrapper(\&Genome::Disk::Volume::percent_used);

    my $super_execute = $self->super_can('_execute_body');
    return $self->$super_execute(@_);
}

sub make_is_mounted_wrapper {
    my($self, $original_sub) = @_;

    return sub {
        my $volume = shift;

        if ($self->is_volume_mounted($volume)) {
            return $volume->$original_sub(@_);
        } else {
            return '<unmounted>';
        }
    };
}

my %cached_is_mounted;
sub is_volume_mounted {
    my($self, $volume) = @_;

    unless (%cached_is_mounted) {
        _populate_cached_is_mounted();
    }

    my $path = $volume->is_remote_volume ? $volume->physical_path : $volume->mount_path;

    if (! exists($cached_is_mounted{$path}) && $self->accurate_size) {
        my $dir = IO::Dir->new($path);
        $dir->read if $dir;
        _populate_cached_is_mounted();
        $cached_is_mounted{$path} //= 0;
    }

    return $cached_is_mounted{$path}
}

sub _populate_cached_is_mounted {
    my $mounts = IO::File->new('/proc/mounts');
    while(my $line = $mounts->getline) {
        my(undef, $mount_path) = split(/\s/, $line);
        $cached_is_mounted{$mount_path} = 1;
    }
}

1;

