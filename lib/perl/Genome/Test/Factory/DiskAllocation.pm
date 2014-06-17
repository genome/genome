package Genome::Test::Factory::DiskAllocation;

use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

require File::Temp;
require File::Path;
use Genome;
require Sub::Install;

my $tmp_volume;
sub tmp_volume {
    return $tmp_volume if $tmp_volume;
    $tmp_volume = Genome::Disk::Volume->__define__(
        mount_path => File::Temp::tempdir(CLEANUP => 1),
        disk_status => 'active',
        can_allocate => 1,
        total_kb => 1073741824,
    ),
    return $tmp_volume;
}

my $archived_volume;
sub archived_volume {
    return $archived_volume if $archived_volume;
    Genome::Disk::Volume->__define__(
        mount_path => Genome::Disk::Volume->archive_volume_prefix,
        disk_status => 'active',
        can_allocate => 1,
        total_kb => 1073741824,
    );
    return $archived_volume;
}

sub generate_obj {
    my ($self, %params) = @_;

    my $owner = delete $params{owner};
    Carp::confess('No owner given!') if not $owner;
    $params{owner_id} = $owner->id;
    $params{owner_class_name} = $owner->class;

    $params{mount_path} //= $self->tmp_volume->mount_path;
    $params{kilobytes_requested} //= 2000;
    $params{disk_group_name} //= 'info_apipe';
    $params{group_subdirectory} //= Genome::Test::Factory::Util::generate_name('group_subdirectory');
    $params{allocation_path} //= Genome::Test::Factory::Util::generate_name('allocation_path');

    my $disk_allocation = Genome::Disk::Allocation->__define__(%params);
    File::Path::mkpath($disk_allocation->absolute_path);

    Sub::Install::reinstall_sub({
            code => sub{ $_[0]->mount_path( $self->archived_volume->mount_path ); return 1; },
            into => $disk_allocation, #'Genome::Disk::Allocation',
            as   => 'archive',
        });

    Sub::Install::reinstall_sub({
            code => sub{ $_[0]->mount_path( $self->tmp_volume()->mount_path ); return 1; },
            into => $disk_allocation, #'Genome::Disk::Allocation',
            as   => 'unarchive',
        });

    return $disk_allocation;
}

1;

