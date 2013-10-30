package Genome::Disk::Group;

use strict;
use warnings;

class Genome::Disk::Group {
    table_name => 'disk.group',
    id_by => [
        dg_id => { is => 'Number', column_name => 'id' },
    ],
    has => [
        disk_group_name => { is => 'Text', column_name => 'name' },
        permissions => { is => 'Number' },
        setgid => { is => 'Number', is_transient => 1 }, #transient during transition from "sticky" column
        subdirectory => { is => 'Text' },
        unix_uid => { is => 'Number' },
        unix_gid => { is => 'Number' },
        user_name => {
            calculate_from => 'unix_uid',
            calculate => q|
                my ($user_name) = $self->_resolve_user_name($unix_uid);
                return $user_name;
            |,
        },
        group_name => {
            calculate_from => 'unix_gid',
            calculate => q| 
                my ($group_name) = $self->_resolve_group_name($unix_gid);
                return $group_name;
            |,
        },
    ],
    has_many_optional => [
        mount_paths => {
            via => 'volumes',
            to => 'mount_path',
        },
        volumes => {
            is => 'Genome::Disk::Volume',
            via => 'assignments',
            to =>  'volume',
        },
        assignments => {
            is => 'Genome::Disk::Assignment',
            reverse_id_by => 'group',
        },
    ],
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'Represents a disk group (eg, info_apipe), which contains any number of disk volumes',
};


my %user_name_cache;
my %group_name_cache;

# memoizing frontends for user_name and group_name
sub _resolve_user_name {
    my($self, $uid) = @_;

    unless (exists $user_name_cache{$uid}) {
        ($user_name_cache{$uid}) = getpwuid($uid);
    }
    return $user_name_cache{$uid};
}

sub _resolve_group_name {
    my($self,$gid) = @_;

    unless (exists $group_name_cache{$gid}) {
        ($group_name_cache{$gid}) = getgrgid($gid);
    }
    return $group_name_cache{$gid};
}


1;
