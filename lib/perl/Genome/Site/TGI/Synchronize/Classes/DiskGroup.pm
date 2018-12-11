package Genome::Site::TGI::Synchronize::Classes::DiskGroup;

use strict;
use warnings;

use Genome;
use Carp qw(confess);
use Memoize qw(memoize);;
use Module::Find qw(findsubmod usesub);

usesub Genome::Disk::Group::Validate;

class Genome::Site::TGI::Synchronize::Classes::DiskGroup {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'DISK_GROUP',
    id_by => [
        id => { is => 'Number', column_name => 'DG_ID' },
    ],
    has => [
        name => { is => 'Text', column_name => 'DISK_GROUP_NAME' },
        permissions => { is => 'Number' },
        setgid => { is => 'Number', is_transient => 1, is_optional => 1 },
        subdirectory => { is => 'Text' },
        unix_uid => { is => 'Number' },
        unix_gid => { is => 'Number' },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

sub genome_class_for_create { return 'Genome::Disk::Group' } 

sub entity_name { return 'disk group'; }

sub properties_to_copy {
    return (qw(
        id
        name
        permissions
        setgid
        subdirectory
        unix_uid
        unix_gid
    ));
}

1;
