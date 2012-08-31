package Genome::ModelGroup::Command::SetDoNotArchive;

use strict;
use warnings;
use Genome;

class Genome::ModelGroup::Command::SetDoNotArchive {
    is => 'Command::V2',
    has => [
        model_groups => {
            is => 'Genome::ModelGroup',
            is_many => 1,
            shell_args_position => 1,
            doc => 'model groups whose members\' last complete build are not to be archived',
        },
    ],
    has_optional => [
        reason => {
            is => 'Text',
            doc => 'reason for marking the builds',
        },
    ],
};

sub help_detail {
    return 'The last complete build of each model in the provided model groups is marked so they can\'t ' .
        'be archived. This process marks all inputs and software results used by the builds, preventing them ' .
        'from ever being archived';
}
sub help_brief { return 'last complete build of each model in the provided group(s) marked so they can\'t be archived' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    for my $group ($self->model_groups) {
        next unless $group;
        for my $build (map { $_->last_complete_build } $group->models) {
            next unless $build;
            for my $allocation ($build->all_allocations) {
                next unless $allocation;
                next unless $allocation->archivable;
                $allocation->archivable(0, $self->reason);
            }
        }
    }
    $self->status_message("All last complete builds of models in provided groups marked, will not be archived");
    return 1;
}

1;

