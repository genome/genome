package Genome::Model::Build::Command::DetermineError;

use strict;
use warnings;

use Genome;

sub handle_failed_from_logs {
    my ($self, $workflow) = @_;

    my @failed_steps = failed_workflow_steps($workflow);
    for my $failed_step (@failed_steps) {
        if ($failed_step->current->can('stderr') and my $error_log = $failed_step->current->stderr) {
            if (-e $error_log and -s $error_log) {

                $self->set_status_from_log_file($error_log);

                # if we can't determine the date from the log file, try the workflow step.
                unless ($self->error_date) {
                    if ($failed_step->end_time) {
                        $self->error_date($failed_step->end_time);
                    } elsif ($failed_step->start_time) {
                        $self->error_date("Sometime After " . $failed_step->start_time);
                    }
                }
                return;
            }
        }
    }
}

sub failed_workflow_steps {
    my $workflow = shift;

    my $failed_steps = [];
    workflow_visitor($workflow, $failed_steps);

    return sort {
        normalize_parent_id($a->parent_instance_id) <=> normalize_parent_id($b->parent_instance_id) ||
        step_time($a) cmp step_time($b)
    } @$failed_steps;
}

sub step_time {
    my $step = shift;
    return $step->end_time if $step->end_time;
    return 's' . $step->start_time if $step->start_time;
    return 'Unknown';
}

sub normalize_parent_id {
    my $parent_id = shift;
    return $parent_id ? -1 : 1;
}

sub workflow_visitor {
    my ($step, $failed_steps) = @_;

    my $status = $step->status;

    if($status eq 'failed' or $status eq 'crashed' or
       $status eq 'running' or $status eq 'running*') {
        push(@{$failed_steps}, $step);
    }

    if ($step->can('related_instances')) {
        for my $sub_step ($step->related_instances) {
            workflow_visitor($sub_step, $failed_steps);
        }
    }
}

1;
