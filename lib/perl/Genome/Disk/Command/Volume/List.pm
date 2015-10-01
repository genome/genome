package Genome::Disk::Command::Volume::List;

use strict;
use warnings;

use Genome;

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
    ],
    doc => 'Lists Genome::Disk::Volume objects',
};

sub execute {
    my $self = shift;

    require Genome::Disk::Group;
    require Genome::Disk::Volume;
    no warnings 'redefine';
    local *Genome::Disk::Group::validate = sub { 1 };  # For just listing, don't validate them
    local *Genome::Disk::Volume::used_kb = make_is_mounted_wrapper(\&Genome::Disk::Volume::used_kb);
    local *Genome::Disk::Volume::percent_used = make_is_mounted_wrapper(\&Genome::Disk::Volume::percent_used);

    my $super_execute = $self->super_can('_execute_body');
    return $self->$super_execute(@_);
}

my %cached_is_mounted;
sub make_is_mounted_wrapper {
    my $original_sub = shift;

    return sub {
        my $volume = shift;

        if ($cached_is_mounted{$volume->id} //= $volume->is_mounted) {
            return $volume->$original_sub(@_);
        } else {
            return '<unmounted>';
        }
    };
}

1;

