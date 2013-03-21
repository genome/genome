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

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

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

1;
