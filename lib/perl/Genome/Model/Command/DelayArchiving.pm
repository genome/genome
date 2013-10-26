package Genome::Model::Command::DelayArchiving;

use strict;
use warnings;
use Genome;

class Genome::Model::Command::DelayArchiving {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            doc => 'models to have archiving of their latest complete build delayed',
        },
    ],
    has_optional => [
        Genome::Disk::Command::Allocation::DelayArchiving->_duration_property_for_commands()
    ],
};

sub help_detail {
    return 'Archiving for the last complete build for each model is pushed into the future. ' .
        'This process marks all allocations those builds use, including instrument data, ' .
        'alignment results, variation detection results, etc';
}
sub help_brief { return 'delays archiving of the last completed build of each model' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    my @builds;
    for my $model ($self->models) {
        push @builds, $model->last_complete_build;
    }
    Genome::Model::Build::Command::DelayArchiving->execute(
        builds => [grep {$_} @builds],
        duration => $self->duration,
    );
    $self->status_message("Last complete build of provided models wont be archived for awhile!");
    return 1;
}

1;

