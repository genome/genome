package Genome::Config::Set;

use strict;
use warnings;

class Genome::Config::Set {
    is => 'Genome::Utility::ObjectWithTimestamps',
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'GENOME_CONFIG_SET',
    id_by => [
        id => {
            is => 'Text',
        },
    ],
    has => [
        allocation_id => {
            is => 'Text',
            column_name => 'allocation_id'
        },
        allocation => {
            is => 'Genome::Disk::Allocation',
            id_by => 'allocation_id',
        },
        path => {
            via => 'allocation',
            to => 'absolute_path',
            is_optional => 1,
        }
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    eval {
        $self->_create_allocation();
    };
    if(my $error = $@) {
        $self->delete();
        die($error);
    }
    return $self;
}

sub _create_allocation {
    my $self = shift;
    unless($self->allocation){
        my $allocation = Genome::Disk::Allocation->create(
            owner_id            => $self->id,
            disk_group_name     => 'info_apipe_ref',
            allocation_path     => 'analysis_configuration/' . $self->id,
            owner_class_name    => 'Genome::Config::Set',
            kilobytes_requested => 25,
        );

        $self->allocation($allocation);
    }
}

1;
