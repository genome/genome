
# review gsanders jlolofie
# note: maybe calculate usage estmate instead of hardcoded value

package Genome::Model::Build::Somatic;
#:adukes this looks fine, there may be some updates required for changes to model inputs and new build protocol, ebelter will be a better judge

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Somatic {
    is => 'Genome::Model::Build',
    has_optional => [
        tumor_build_links                  => { is => 'Genome::Model::Build::Link', reverse_as => 'to_build', where => [ role => 'tumor'], is_many => 1,
                                               doc => 'The bridge table entry for the links to tumor builds (should only be one)' },
        tumor_build                     => { is => 'Genome::Model::Build', via => 'tumor_build_links', to => 'from_build',
                                               doc => 'The tumor build with which this build is associated' },
        normal_build_links                  => { is => 'Genome::Model::Build::Link', reverse_as => 'to_build', where => [ role => 'normal'], is_many => 1,
                                               doc => 'The bridge table entry for the links to normal builds (should only be one)' },
        normal_build                     => { is => 'Genome::Model::Build', via => 'normal_build_links', to => 'from_build',
                                               doc => 'The tumor build with which this build is associated' },
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
    return $self->build_id . ' Somatic Pipeline';
}

# Returns the newest somatic workflow instance associated with this build
# Note: Only somatic builds launched since this code was added will have workflows associated in a queryable manner
sub newest_somatic_workflow_instance {
    my $self = shift;

    my $build_workflow = $self->newest_workflow_instance;
    return unless $build_workflow;

    if($build_workflow->name =~ 'Somatic Pipeline') {
        #Newer builds run the pipeline's workflow directly
        return $build_workflow;
    }

    #Older builds had many layers of indirection that eventually lead to a workflow with this name
    my @sorted = sort {
        $b->id <=> $a->id
    } Workflow::Operation::Instance->get(
        name => 'Somatic Pipeline Build ' . $self->build_id
    );

    unless (@sorted) {
        $self->warning_message("No somatic workflow instances found for build " . $self->id);
        return;
    }

    return $sorted[0];
}

1;
