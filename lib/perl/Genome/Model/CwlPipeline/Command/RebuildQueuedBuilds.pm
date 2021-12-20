package Genome::Model::CwlPipeline::Command::RebuildQueuedBuilds;

use strict;
use warnings;

use Genome;

class Genome::Model::CwlPipeline::Command::RebuildQueuedBuilds {
    is => 'Command::V2',
    roles => ['Genome::Model::Command::Submittable'],
    doc => 'process the rebuild queue list for the current user',
    has => [
        limit => {
            doc => 'only restart up to this many builds',
            is => 'Integer',
            default => 10,
        },
        submit_jobs => {
            doc => 'if set, will use `genome analysis-project possible-volumes` to submit jobs for the build starts',
            is => 'Boolean',
            default => 0,
        },
    ],
};

sub help_detail {
    return <<EOS;
Looks for queued restart requests for the current user and performs the restarts.
EOS
}

sub execute {
    my $self = shift;

    my @builds_to_restart = Genome::Model::Build::CwlPipeline->get(
        run_by => Genome::Sys->username,
        action_requested => 'restart',
        status => 'Failed',
        -order => ['created_at'],
        -limit => $self->limit,
    );

    my %builds_by_anp;

    for my $build (@builds_to_restart) {
        my $anp = $build->model->analysis_project;
        $builds_by_anp{$anp->id} //= [];
        push @{ $builds_by_anp{$anp->id} }, $build;
    }

    my $restart_count = 0;
    for my $anp_id (keys %builds_by_anp) {
        my $anp = Genome::Config::AnalysisProject->get($anp_id);

        my $builds = $builds_by_anp{$anp_id};

        if ($self->submit_jobs) {
            $self->_submit_jobs(
                $anp,
                [qw(genome model cwl-pipeline restart), map $_->id, @$builds],
            );
        } else {
            my $cmd = Genome::Model::CwlPipeline::Command::Restart->create(
                builds => $builds,
            );
            $cmd->execute();
        }

        $restart_count += scalar @$builds;
    }

    my $remaining_count = Genome::Model::Build::CwlPipeline->define_set(
        run_by => Genome::Sys->username,
        action_requested => 'restart',
        status => 'Failed',
    )->count;
    $self->status_message('Restarted %s builds; %s left in queue', $restart_count, $remaining_count);

    return 1;
}

1;
