package Genome::InstrumentData::Sanger;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Sanger {
    is => 'Genome::InstrumentData',
    has_constant => [
        sequencing_platform => { value => 'sanger' },
        read_count => { calculate => q|return 96;|,},
    ],
    has_optional => [
        research_project => {
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'research_project' ],
            is_mutable => 1,
            default_value => 'unknown',
        },
        disk_allocation => {
            is => 'Genome::Disk::Allocation',
            calculate_from => [ 'subclass_name', 'id' ],
            calculate => q{ return Genome::Disk::Allocation->get(owner_id => $id, owner_class_name => $subclass_name); },
        },
    ],
};

sub full_path {
    my $self = shift;

    my $disk_allocation = $self->disk_allocation;
    return if not $disk_allocation;

    return $disk_allocation->absolute_path;
}

sub dump_to_file_system {
    my $self = shift;

    my $rv = eval { 
        Genome::Sys->shellcmd(
            cmd => 'gmt lims import-sanger-runs '.$self->id,
        );
    };
    return 1 if $rv;

    $self->error_message('Failed to dump reads to file system for run '.$self->id);
    return;
}

1;

