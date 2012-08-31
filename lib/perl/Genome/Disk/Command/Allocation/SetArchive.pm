package Genome::Disk::Command::Allocation::SetArchive;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::SetArchive {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Allocations that are to be set as archivable',
        },
    ],
    has_optional => [
        reason => {
            is => 'Text',
            doc => 'reason for wanting to set these allocations as archivable',
        },
    ],
};

sub help_detail { return 'given allocations are set as archivable' };
sub help_brief { return help_detail() };
sub help_synopsis { return help_detail() };

sub execute {
    my $self = shift;
    for my $allocation ($self->allocations) {
        next if $allocation->archivable;
        $allocation->archivable(1, $self->reason);
    }
    $self->status_message("Successfully set allocations as archivable");
    return 1;
}

1;

