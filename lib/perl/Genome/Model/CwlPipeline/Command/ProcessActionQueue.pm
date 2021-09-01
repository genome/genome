package Genome::Model::CwlPipeline::Command::ProcessActionQueue;

use strict;
use warnings;

use Genome;

class Genome::Model::CwlPipeline::Command::ProcessActionQueue {
    is => 'Command::V2',
    roles => ['Genome::Model::Command::Submittable'],
    has => [
        max_count => {
            is => 'Number',
            doc => 'Process at most this many queued builds',
            default => 10,
        },
        submit_jobs => {
            is => 'Boolean',
            default => 0,
            doc => 'if set, will use `genome analysis-project possible-volumes` to submit jobs for the actions',
        },
    ],
};

sub help_detail {
    return <<EOHELP
Takes the actions requested for CWLPipeline builds.
EOHELP
}

sub execute {
    my $self = shift;

    my @to_process = Genome::Model::Build::CwlPipeline->get(
        action_requested => ['stop', 'abandon'],
        -limit => $self->max_count,
    );

    my %builds_by_anp;

    for my $build (@to_process) {
        my $anp = $build->model->analysis_project;
        my $action = $build->action_requested;

        $builds_by_anp{$anp->id}{$action} //= [];
        push @{ $builds_by_anp{$anp->id}{$action} }, $build;
    }

    for my $anp_id (keys %builds_by_anp) {
        my $anp = Genome::Config::AnalysisProject->get($anp_id);

        for my $action (keys %{ $builds_by_anp{$anp_id} }) {
            my $builds = $builds_by_anp{$anp_id}{$action};

            if ($self->submit_jobs) {
                $self->_submit_jobs(
                    $anp,
                    [qw(genome model build), $action, map $_->id, @$builds],
                );
            } else {
                for my $build (@$builds) {
                    $build->$action;
                }
            }

            for my $build (@$builds) {
                $build->action_requested(0);
            }
        }
    }

    return 1;
}

1;
