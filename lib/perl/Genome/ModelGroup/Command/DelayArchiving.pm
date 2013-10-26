package Genome::ModelGroup::Command::DelayArchiving;

use strict;
use warnings;
use Genome;

class Genome::ModelGroup::Command::DelayArchiving {
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
        Genome::Disk::Command::Allocation::DelayArchiving->_duration_property_for_commands()
    ],
};

sub help_detail {
    return 'Archiving of the last complete build of each model in the provided model groups is delayed ' .
        'This process marks all inputs and software results used by the builds, preventing them ' .
        'from being archived for awhile';
}
sub help_brief { return 'prevents the last completed build from models in the model group from being archived' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    my @models;
    for my $group ($self->model_groups) {
        push @models, $group->models
    }
    Genome::Model::Command::DelayArchiving->execute(
        duration => $self->duration,
        models => \@models,
    );
    $self->status_message("All last complete builds of models in provided groups marked, will not be archived");
    return 1;
}

1;

