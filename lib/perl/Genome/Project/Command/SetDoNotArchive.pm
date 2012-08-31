package Genome::Project::Command::SetDoNotArchive;

use strict;
use warnings;
use Genome;

class Genome::Project::Command::SetDoNotArchive {
    is => 'Command::V2',
    has => [
        projects => {
            is => 'Genome::Project',
            is_many => 1,
            shell_args_position => 1,
            doc => 'projects that are not to be archived',
        },
    ],
    has_optional => [
        reason => {
            is => 'Text',
            doc => 'reason provided projects are not to be archived',
        },
    ],
};

sub help_detail {
    return 'Any allocations associated with the given projects  will have their archivable ' .
        'flag set to false, preventing them from being archived.';
}
sub help_brief { return 'marks allocations associated with given projects as unable to be archived' };
sub help_synopsis { return help_brief . "\n" };

sub execute {
    my $self = shift;
    for my $project ($self->projects) {
        next unless $project;

        # For now, it is safe to assume that the only thing in a project that matters are the models.
        # Anything else is ignored.

        my @models = grep { $_->isa('Genome::Model') } $project->entities;
        for my $model (@models) {
            my $build = $model->last_complete_build;
            next unless $build;
            for my $allocation ($build->all_allocations) {
                next unless $allocation;
                next unless $allocation->archivable;
                $allocation->archivable(0, $self->reason);
            }
        }
    }

    $self->status_message("Last complete build of all models in provided projects marked as do-not-archive");
    return 1;
}

1;

