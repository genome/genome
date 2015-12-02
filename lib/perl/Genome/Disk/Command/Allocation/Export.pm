package Genome::Disk::Command::Allocation::Export;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::Export {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Allocations to export',
        },
        output_dir => {
            is => 'Text',
            doc => 'The output directory to copy allocation contents to.',
        },
    ],
    doc => 'copy allocations from allocated disk to another network filesystem path.',
};

sub help_detail {
    return <<EOS
Copies allocations from allocated disk to a specified output directory. 

EOS
}

sub help_brief {
    return 'copy allocations from allocated disk to a specified output directory';
}

sub execute {
    my $self = shift;

    my $existing_allocation = Genome::Disk::Allocation->get_allocation_for_path($self->output_dir);
    if ( $existing_allocation && ($existing_allocation->absolute_path eq $self->output_dir) ) {
        die('Can not export allocation to allocated disk directory: '. $self->output_dir);
    }

    for my $allocation ($self->allocations) {
        $allocation->copy( output_dir => File::Spec->join($self->output_dir, $allocation->id) );
    }

    $self->status_message("Successfully exported allocations!");
    return 1;
}

1;

