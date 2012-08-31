package Genome::Task::Command::Run;

use warnings;
use strict;
use Genome;
use JSON::XS;
use File::Path;

class Genome::Task::Command::Run {
    is => 'Command::V2',
    has => [
        task => {
            is => 'Genome::Task',
            doc => 'Task to run.  Resolved from command line via text string',
            shell_args_position => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    
    my $class_object;
    eval {$class_object = $self->task->command_class->class};
    if (!$class_object) {
        $self->error_message(sprintf("Failed to execute, %s couldn't be loaded", $self->task->command_class));
        return $self->handle_failure;
    }

    my $cmd_object = $self->task->command_object;
    my $task_id = $self->task->id;
    if (!$cmd_object) {
        $self->error_message(sprintf("Failed to execute, %s couldn't be instantiated", $self->task->command_class));
        return $self->handle_failure;
    }

    my %attrs_to_update = (status => 'running', time_started => $UR::Context::current->now);

    my $transaction = UR::Context::Transaction->begin;
    $self->task->out_of_band_attribute_update(%attrs_to_update);
#    $self->task->unload();
#    $self->task(UR::Context->current->reload('Genome::Task', id=>$task_id));


    my $result;
    eval {
        $result = $cmd_object->execute;
    };
    
    if ($@ || !$result) {
        $self->error_message("COMMAND FAILURE:  $@ -- " . $cmd_object->error_message);
        $transaction->rollback;
        $self->task->out_of_band_attribute_update(status=>'failed');
#        $self->task->unload();
#        $self->task(UR::Context->current->reload('Genome::Task', id=>$task_id));
    } else {
        $transaction->commit;
        $self->task->status("succeeded");
    }

    $self->task->time_finished($UR::Context::current->now);

    return $result;
}

1;
