package Genome::Model::Command::Services::Build::BmodWorkflowPriority;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::Command::Services::Build::BmodWorkflowPriority { 
    is => 'Genome::Command::Base',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            is_many => 1,
            doc => 'Builds to modify the workflow jobs for.',
        },
        priority => {
            is => 'Integer',
            doc => 'LSF priority for job, 0..500',
            default => '500',
        },
    ],
};

sub execute {
    my $self = shift;
    my @builds = $self->builds;
    my $priority = $self->priority;

    for my $build (@builds) {
        my $job_id = $self->job_id_for_build($build);
        unless ($job_id) {
            $self->status_message("No LSF job for build (" . $build->id . ")");
            next;
        }

        my $host_group = $self->parse_host_group_of_lsf_job($job_id);

        $self->status_message("Setting user job priority to $priority for LSF job $job_id.");
        system("bmod -sp $priority -m $host_group $job_id");
    }

    return 1;
}

sub job_id_for_build {
    my $self = shift;
    my $build = shift;

    my $job_id;
    my $workflow = $build->newest_workflow_instance;
    $job_id = $workflow->current->dispatch_identifier unless ($job_id || !$workflow);
    $job_id = $build->the_master_event->lsf_job_id unless ($job_id);

    return $job_id;
}

sub parse_host_group_of_lsf_job {
    my $self = shift;
    my $job_id = shift;

    chomp(my $queue = qx(bjobs $job_id | grep ^[0-9] | awk '{print \$4}'));
    my $bqueus_cmd_output = qx(bqueues -l $queue | grep ^HOSTS:);
    my ($host_group) = $bqueus_cmd_output =~ /^HOSTS:\s+(\w+)\//;

    return $host_group;
}

1;
