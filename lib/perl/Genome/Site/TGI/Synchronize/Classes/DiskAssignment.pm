package Genome::Site::TGI::Synchronize::Classes::DiskAssignment;
use Genome;
use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::DiskAssignment {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
    table_name => 'DISK_VOLUME_GROUP',
    id_by => [
        group_id => {
            is => 'Number',
            column_name => 'DG_ID',
        },
        volume_id => {
            is => 'Number',
            column_name => 'DV_ID',
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

sub genome_class_for_create { return 'Genome::Disk::Assignment' }

sub entity_name { return 'disk assignment'; }

sub properties_to_copy {
    return (qw(
        group_id
        volume_id
    ));
}

1;
