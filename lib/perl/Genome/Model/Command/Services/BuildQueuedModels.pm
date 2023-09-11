package Genome::Model::Command::Services::BuildQueuedModels;

use strict;
use warnings;

use Genome;

use Try::Tiny qw(try catch);

class Genome::Model::Command::Services::BuildQueuedModels {
    is => 'Command::V2',
    roles => ['Genome::Model::Command::Submittable'],
    doc => 'process the build queue list for the current user',
    has => [
        scheduled_max => {
            doc => 'maximum number of scheduled builds',
            is => 'Integer',
            default => 50,
        },
        per_analysis_project_scheduled_max => {
            doc => 'maximum number of scheduled builds for a single analysis project',
            is => 'Integer',
            default => 10,
        },
        running_max => {
            doc => 'maximum number of running builds',
            is => 'Integer',
            default => 750,
        },
        start_limit => {
            doc => 'only start up to this many builds even if more are in the queue and would fit',
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
Looks for queued builds for the current user to start and starts them if there are available slots.
EOS
}

sub execute {
    my $self = shift;

    my $running_count = $self->_current_build_count('Running');
    if ($running_count >= $self->running_max) {
        $self->status_message(
            'Running build count %s meets or exceeds maximum %s.',
            $running_count, $self->running_max,
        );

        $self->_final_status(0);
        return 1;
    }

    my $scheduled_count = $self->_current_build_count('Scheduled');
    if ($scheduled_count >= $self->scheduled_max) {
        $self->status_message(
            'Scheduled build count %s meets or exceeds maximum %s.',
            $scheduled_count, $self->scheduled_max,
        );

        $self->_final_status(0);
        return 1;
    }

    my $iterator = Genome::Model->create_iterator(
        run_as => Genome::Sys->username,
        build_requested => 1,
        -order => ['creation_date'],
    );

    my %full_anps;
    my $num_started = 0;
    MODEL: while ($num_started < $self->start_limit and $scheduled_count < $self->scheduled_max and my $model = $iterator->next) {
        my $anp = $model->analysis_project;

        next MODEL if $full_anps{$anp->id}; #already found this AnP was full

        my $per_anp_scheduled_count = $self->_current_build_count('Scheduled', $anp->id);
        if ($per_anp_scheduled_count >= $self->per_analysis_project_scheduled_max) {
            $self->status_message(
                'Scheduled build count %s for AnP %s meets or exceeds maximum %s',
                $per_anp_scheduled_count, $anp->id, $self->per_analysis_project_scheduled_max,
            );
            $full_anps{$anp->id} = 1;
            next MODEL;
        }

        if ($self->submit_jobs) {
            try {
                $self->_submit_jobs(
                    $anp,
                    [qw(genome model build start --force --unstartable-ok), $model->id],
                );
            } catch {
                my $error = $_;
                $self->error_message('Failed to bsub build start for %s: %s', $model->id, $error);
            };
        } else {
            my $cmd = Genome::Model::Build::Command::Start->create(
                models => [$model],
                force => 1,
                unstartable_ok => 1,
            );
            $cmd->execute();
        }

        $num_started++;
        $scheduled_count++;
    }

    undef $iterator;

    $self->_final_status($num_started);

    return 1;
}

sub _current_build_count {
    my $self = shift;
    my $status = shift;
    my $anp_id = shift;

    my @query = (run_by => Genome::Sys->username, status => $status);
    if ($anp_id) {
        push @query,
            'model.analysis_project.id' => $anp_id;
    }

    return Genome::Model::Build->define_set(@query)->count;
}

sub _final_status {
    my $self = shift;
    my $num_started = shift;

    my $remaining = Genome::Model->define_set(build_requested => 1, run_as => Genome::Sys->username)->count;

    $self->status_message('Started %s builds; %s left in queue.', $num_started, $remaining);
}

1;

