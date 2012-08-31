package Genome::Model::Build::Command::SetDoNotArchive;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Command::SetDoNotArchive {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'builds are marked so they cannot be archived',
        },
    ],
    has_optional => [
        reason => {
            is => 'Text',
            doc => 'reason for marking these builds'
        },
    ],
};

sub help_detail {
    return 'All allocations that are used by this build are marked such that they can\'t be archived. ' .
        'This includes allocations the build itself owns, all allocations that its ' .
        'inputs own, and allocations owned by software results this build uses.'
}
sub help_brief { return 'allocations used by provided builds are marked so they can\'t be archived' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    for my $build ($self->builds) {
        for my $allocation ($build->all_allocations) {
            next unless $allocation->archivable;
            $allocation->archivable(0, $self->reason);
        }
    }
    $self->status_message("Provided builds now can't be archived");
    return 1;
}

1;

