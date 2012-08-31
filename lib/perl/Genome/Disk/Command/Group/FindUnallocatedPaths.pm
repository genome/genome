package Genome::Disk::Command::Group::FindUnallocatedPaths;

use warnings;
use strict;

use Genome;

class Genome::Disk::Command::Group::FindUnallocatedPaths{
    is => 'Command::V2',
    has_input => [
        group => {
            is => 'Genome::Disk::Group',
            doc => 'Identifier for group on which to find unallocated paths.',
        },
    ],
    has_optional => [
        _unallocated_paths => {
            is => 'Text',
            is_many => 1,
            is_output => 1,
        },
    ],
    doc => 'Finds unallocated paths in the specified disk group',
};

sub help_detail {
    return 'Scans volumes in the provided group for paths that are not contained in an allocation';
}

sub execute{
    my $self = shift;
    for my $volume ($self->group->volumes) {
        #we only care for active drives
        unless($volume->disk_status eq 'active') {
            next;
        }
        #and drives whene we might allocate files
        unless($volume->can_allocate){
            next;
        }
        my $cmd = Genome::Disk::Command::Volume::FindUnallocatedPaths->create(
            volume => $volume,
        );
        $cmd->execute;
        my @paths = $self->_unallocated_paths(), $cmd->_unallocated_paths();
        $self->_unallocated_paths(\@paths);
    }
}
