package Genome::Disk::Assignment;
use Genome;
use strict;
use warnings;

class Genome::Disk::Assignment {
    table_name => 'disk.volume_group_bridge',
    id_by => [
        group_id => {
            is => 'Number',
            doc => 'disk group ID',
        },
        volume_id => {
            is => 'Number',
            doc => 'disk volume ID'
        },
    ],
    has => [
        group => {
            is => 'Genome::Disk::Group',
            id_by => 'group_id',
        },
        disk_group_name => { via => 'group' },
        user_name => { via => 'group' },
        group_name => { via => 'group' },
        subdirectory => { via => 'group' },
        volume => {
            is => 'Genome::Disk::Volume',
            id_by => 'volume_id',
        },
        mount_path => { via => 'volume' },
        total_kb   => { via => 'volume' },
        cached_unallocated_kb => { via => 'volume' },
        unallocated_kb => { via => 'volume' },
        percent_allocated => { via => 'volume' },
        percent_used => { via => 'volume' },
        absolute_path => {
            calculate_from => ['mount_path','subdirectory'],
            calculate => q| return $mount_path .'/'. $subdirectory; |,
        },
    ],
    data_source => 'Genome::DataSource::GMSchema',
};

1;
