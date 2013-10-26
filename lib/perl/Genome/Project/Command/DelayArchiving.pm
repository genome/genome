package Genome::Project::Command::DelayArchiving;

use strict;
use warnings;
use Genome;

class Genome::Project::Command::DelayArchiving {
    is => 'Command::V2',
    has => [
        projects => {
            is => 'Genome::Project',
            is_many => 1,
            shell_args_position => 1,
            doc => 'projects to be kept longer',
        },
    ],
    has_optional => [
        Genome::Disk::Command::Allocation::DelayArchiving->_duration_property_for_commands()
    ],
};

sub help_detail {
    return 'Any allocations associated with the given projects will have their archiving date pushed back' .
        'by the duration specified, preventing them from being archived.';
}
sub help_brief { return 'delays archiving of allocations associated with this project' };
sub help_synopsis { return help_brief . "\n" };

sub execute {
    my $self = shift;
    for my $project ($self->projects) {
        next unless $project;

        # For now, it is safe to assume that the only thing in a project that matters are the models.
        # Anything else is ignored.
        my @models = grep { $_->isa('Genome::Model') } $project->entities;
        Genome::Model::Command::DelayArchiving->execute(
            duration => $self->duration,
            models => \@models,
        );
    }

    $self->status_message("Archiving of last complete build of all models in provided projects delayed");
    return 1;
}

1;

