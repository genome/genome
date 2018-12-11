package Genome::Site::TGI::Synchronize::Classes::DiskVolume;

use strict;
use warnings;

use Genome;
use Carp;

use Data::Dumper;
use Filesys::Df qw();
use List::Util qw(max);
use Scope::Guard;

class Genome::Site::TGI::Synchronize::Classes::DiskVolume {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'DISK_VOLUME',
    id_by => [
        id => {is => 'Number', column_name => 'DV_ID'},
    ],
    has => [
        hostname => { is => 'Text' },
        physical_path => { is => 'Text' },
        mount_path => { is => 'Text' },
        disk_status => {
            is => 'Text',
        },
        can_allocate => {
            is => 'Number',
        },
        total_kb => { is => 'Number' },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

sub genome_class_for_create { return 'Genome::Disk::Volume' }

sub entity_name { return 'disk volume'; }

sub properties_to_copy {
    return (qw(
        id
        hostname
        physical_path
        mount_path
        disk_status
        can_allocate
        total_kb
    ));
}

sub properties_to_keep_updated {
    return (qw(
        total_kb
    ));
}

1;
