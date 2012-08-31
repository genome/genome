package Genome::Model::Build::SomaticCapture;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::SomaticCapture {
    is => 'Genome::Model::Build',
    has_optional => [
        tumor_build_links                  => { is => 'Genome::Model::Build::Link', reverse_as => 'to_build', where => [ role => 'tumor'], is_many => 1,
                                               doc => 'The bridge table entry for the links to tumor builds (should only be one)' },
        tumor_build                     => { is => 'Genome::Model::Build', via => 'tumor_build_links', to => 'from_build', 
                                               doc => 'The tumor build with which this build is associated' },
        tumor_model                     => { is => 'Genome::Model', via => 'tumor_build', to => 'model' },
                                           normal_build_links                  => { is => 'Genome::Model::Build::Link', reverse_as => 'to_build', where => [ role => 'normal'], is_many => 1,
                                               doc => 'The bridge table entry for the links to normal builds (should only be one)' },
        normal_build                     => { is => 'Genome::Model::Build', via => 'normal_build_links', to => 'from_build', 
                                               doc => 'The tumor build with which this build is associated' },
        normal_model                     => { is => 'Genome::Model', via => 'normal_build', to => 'model' },
    ],
};


sub workflow_instances {
    my $self = shift;
    my @instances = Workflow::Operation::Instance->get(
        name => $self->workflow_name
    );

    #older builds used a wrapper workflow
    unless(scalar @instances) {
        return $self->SUPER::workflow_instances;
    }

    return @instances;
}

sub workflow_name {
    my $self = shift;

    return $self->build_id . ' Somatic Capture Pipeline';
}

# Returns the newest somatic workflow instance associated with this build
# Note: Only somatic builds launched since this code was added will have workflows associated in a queryable manner
sub newest_somatic_workflow_instance {
    my $self = shift;

    my $build_workflow = $self->newest_workflow_instance;
    return unless $build_workflow;
    
    if($build_workflow->name =~ 'Somatic Capture') {
        #Newer builds run the pipeline's workflow directly
        return $build_workflow;
    }

    #Older builds have a wrapper workflow with a stage that executes the real workflow
    my $somatic_capture = Workflow::Operation::Instance->get(parent_instance_id => $build_workflow->id, 'name LIKE' => '% somaticcapture');
    return unless $somatic_capture;

    my $run_workflow = Workflow::Operation::Instance->get(parent_instance_id => $somatic_capture->id, 'name LIKE' => 'run-workflow %');
    return unless $run_workflow;
    
    my @run_workflow_executions = Workflow::Operation::InstanceExecution->get(instance_id => $run_workflow->id);
    return unless scalar(@run_workflow_executions);
    my @somatic_workflows = Workflow::Operation::Instance->get(parent_execution_id => [map($_->id,@run_workflow_executions)]);

    my @sorted = sort {
        $b->id <=> $a->id
    } @somatic_workflows;

    return unless (scalar @sorted);
    
    return $sorted[0];
}

# Returns a hash ref with all of the inputs of the newest somatic workflow instance
sub somatic_workflow_inputs {
    my $self = shift;

    # TODO Switched to doing a direct database query to find inputs, since if we go through the object layer, workflows with steps which have at some point changed class paths 
    # will crash, with no good solution. The best solution is probably not to query the workflow at all, and instead log it elsewherorkflow_instance_namee
    my $ds = $UR::Context::current->resolve_data_sources_for_class_meta_and_rule(Workflow::Operation::Instance->__meta__);
    my $dbh = $ds->get_default_dbh;
    $dbh->{LongReadLen} = 1024*1024;

    my $workflow_instance = $self->newest_somatic_workflow_instance;
    return unless( $workflow_instance);

    my $input_stored = $dbh->selectrow_arrayref("SELECT input_stored FROM workflow_instance WHERE workflow_instance_id = ?", {}, $workflow_instance->id)->[0];

    unless ($input_stored) {
        $self->error_message("Could not find a workflow instance associated with this build for workflow: " . $workflow_instance->name);
        die;
    }

    my $input = Storable::thaw($input_stored);
    unless ($input) {
        $self->error_message("Could not thaw input hash for workflow instance named: " . $workflow_instance->name);
        die;
    }

    # returns hashref of workflow params like { input => value }
    return $input;  
}

# Input: the name of the somatic workflow input you'd like to know
# Returns: value of one input of the latest somatic workflow instance.
# TODO this will break if the build allocations have moved...so... if we ask the workflow for file locations perhaps we should strip off the path and instead use the build's current data_dir
# we could check to see if it changed, and if it did warn and return
sub somatic_workflow_input {
    my $self = shift;
    my $input_name = shift;

    my $input = $self->somatic_workflow_inputs;

    if ($input) {
        unless (exists $input->{$input_name}) {
            my @valid_inputs = sort(keys %$input);
            my $inputs_string = join(", ", @valid_inputs);
            $self->error_message("Input $input_name does not exist. Valid inputs to query for this build are: \n$inputs_string");
            return;
        }

        unless (defined $input->{$input_name}) {
            $self->error_message("Input $input_name exists, but is not defined for this build. Something may have gone wrong with the build.");
            return;
        }

        return $input->{$input_name};
    } else {
        $self->error_message("There was no workflow instance found for this build.");
        return;
    }

}

1;
