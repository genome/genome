package Genome::Disk::Command::Allocation::PathsToAllocations;

use strict;
use warnings;

use Genome;
use File::Basename qw(basename);

class Genome::Disk::Command::Allocation::PathsToAllocations {
    is => 'Command::V2',
    has => [
        absolute_paths => {
            is => 'DirectoryPath',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Absolute path for which an allocation needs to be found',
        },
        filenames => {
            is => 'Boolean',
            default => 0,
            doc => 'Set to true if you would like the output paths to contain the original file (if provided) rather than just the directory'
        },
    ],
};

sub execute {
    my $self = shift;

    my %allocations;
    for my $given_path ($self->absolute_paths) {
        my $allocation = Genome::Disk::Allocation->get_allocation_for_path($given_path);
        if ($allocation) {
            $allocations{$given_path} = $allocation;
        } else {
            $self->warning_message("Found no allocation for path $given_path");
        }
    }

    for my $given_path (keys %allocations) {
        my $allocation = $allocations{$given_path};
        my $allocation_path = $allocation->absolute_path;
        if ($self->filenames) {
            my $file = basename($given_path);
            if (defined $file) {
                $allocation_path = File::Spec->join($allocation_path, $file);
            }
        }
        print join("\t", $allocation->id, $allocation_path) . "\n";
    }
    return 1;
}

1;

