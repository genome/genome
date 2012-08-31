package Genome::Model::Build::Command::ShowWorkflow;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::Command::ShowWorkflow {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'build whose workflow should be displayed',
        },
    ],
};

sub execute {
    my $self = shift;
    my $build = $self->build;

    my $workflow = $build->newest_workflow_instance;
    my $workflow_id = $workflow->id;

    my $output = `workflow show $workflow_id`;
    $self->status_message($output);
    return 1;
}

1;

