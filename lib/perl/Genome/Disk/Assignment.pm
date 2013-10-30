package Genome::Disk::Assignment;
use Genome;
use strict;
use warnings;

class Genome::Disk::Assignment {
    table_name => 'disk.volume_group_bridge',
    id_by => [
        dg_id => {
            is => 'Number',
            column_name => 'group_id',
            doc => 'disk group ID',
        },
        dv_id => {
            is => 'Number',
            column_name => 'volume_id',
            doc => 'disk volume ID'
        },
    ],
    has => [
        group => {
            is => 'Genome::Disk::Group',
            id_by => 'dg_id',
        },
        disk_group_name => { via => 'group' },
        user_name => { via => 'group' },
        group_name => { via => 'group' },
        subdirectory => { via => 'group' },
        volume => {
            is => 'Genome::Disk::Volume',
            id_by => 'dv_id',
        },
        mount_path => { via => 'volume' },
        total_kb   => { via => 'volume' },
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

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    $self->initialize_base_directory;
    return $self;
};

sub initialize_base_directory {
    my $self = shift;
    my $group = $self->group;
    my $volume = $self->volume;
    $volume->mount; 
    
    my $subdirectory = $volume->mount_path . '/' . $group->subdirectory;
    unless (-e $subdirectory) {
        Genome::Sys->create_directory($subdirectory);
    }

}

1;
