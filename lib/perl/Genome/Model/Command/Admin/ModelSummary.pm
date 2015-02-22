package Genome::Model::Command::Admin::ModelSummary;

use strict;
use warnings;

use Genome;

use Try::Tiny qw(try catch);

class Genome::Model::Command::Admin::ModelSummary {
    is => 'Command::V2',
    doc => 'Tool for the Cron Tzar to review builds, e.g. missing builds, failed builds, etc.',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
        },
        hide_statuses => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'Hide build details for statuses listed.',
        },
        hide_no_action_needed => {
            is => 'Boolean',
            default => 0,
            doc => 'Hide models that do not require any action.',
        },
        auto => {
            is => 'Boolean',
            default => 0,
            doc => 'auto act based on recommended actions',
        },
        auto_batch_size  => {
            is => 'Integer',
            default => 100,
            doc => 'When "--auto" is true, commit() to the database after this many Models have been identified for cleanup or restart.  0 means wait until the end to commit everything',
        }
    ],
};


sub execute {
    my $self = shift;

    my @models = $self->models;
    my @hide_statuses = $self->hide_statuses;

    #preload data
    my %failed_build_params = (
        'status in' => ['Unknown', 'Unstartable', 'Failed'],
    );
    my @builds = Genome::Model::Build->get(
        model_id => [map($_->id, @models)],
        %failed_build_params,
    );

    my $synchronous = ($self->auto and $self->auto_batch_size);  # whether we should start builds as we go or wait until the end

    # Header for the report produced at the end of the loop
    $self->print_message(join("\t", qw(model_id action latest_build_status first_nondone_step latest_build_rev model_name pp_name fail_count)));
    $self->print_message(join("\t", qw(-------- ------ ------------------- ------------------ ---------------- ---------- ------- ----------)));

    my %classes_to_unload;
    my $build_requested_count = 0;
    my %cleanup_rv;
    my $change_count = 0;

    my $auto_batch_size_txn = UR::Context::Transaction->begin;
    my $commit = sub {
        unless ($auto_batch_size_txn && $auto_batch_size_txn->isa('UR::Context::Transaction')) {
            die "Not in a transaction! Something went wrong.";
        }

        if($auto_batch_size_txn->commit) {
            $self->debug_message('Committing...');
            unless (UR::Context->commit) {
                die "Commit failed! Bailing out.";
            }

            $self->debug_message('Unloading...');
            # Unload to free up memory.
            for my $class (keys %classes_to_unload) {
                $class->unload;
            }

            $change_count = 0;
            $auto_batch_size_txn = UR::Context::Transaction->begin;
        } else {
            $auto_batch_size_txn->rollback;
        }
    };

    for my $model (@models) {
        my $per_model_txn = UR::Context::Transaction->begin();

        my $build_iterator = $model->build_iterator(
            'status not like' => 'Abandoned',
            '-order_by' => '-created_at',
        );
        my $latest_build        = $build_iterator->next;
        my $latest_build_status = ($latest_build ? $latest_build->status : '-');
        $latest_build_status = 'Requested' if $model->build_requested;


        my @failed_builds = $model->builds(%failed_build_params);
        my $fail_count   = scalar @failed_builds;
        my $model_id     = ($model ? $model->id                       : '-');
        my $model_name   = ($model ? $model->name                     : '-');
        my $pp_name      = ($model ? $model->processing_profile->name : '-');

        my $first_nondone_step;
        if ($latest_build) {
           $first_nondone_step = find_first_nondone_step($latest_build);
        }
        $first_nondone_step ||= '-';

        my $track_change = sub {
            $change_count++;
            $classes_to_unload{$latest_build->class}++ if $latest_build;
        };

        my $request_build = sub {
            $model->build_requested(1, 'ModelSummary recommended rebuild.');
            $build_requested_count++;
        };

        $first_nondone_step =~ s/^\d+\s+//;
        $first_nondone_step =~ s/\s+\d+$//;

        my $latest_build_revision = $latest_build->software_revision if $latest_build;
        $latest_build_revision ||= '-';
        ($latest_build_revision) =~ /\/(genome-[^\/])/ if $latest_build_revision =~ /\/genome-[^\/]/;

        $model_name =~ s/\.?$pp_name\.?/.../;

        my $action;
        if (!$latest_build && $model->build_requested){
            $action = 'none';
        }
        elsif (!$latest_build) {
            $action = 'build-needed';
            if ($self->auto) {
                $request_build->();
                $track_change->();
            }
        }
        elsif ($latest_build_status eq 'Scheduled' || $latest_build_status eq 'Running' || $latest_build_status eq 'Requested') {
            $action = 'none';
        }
        elsif ($latest_build && $latest_build_status eq 'Succeeded' && $fail_count) {
            $action = 'cleanup';
            if ($self->auto) {
                my $cleanup_succeeded = Genome::Model::Command::Admin::CleanupSucceeded->create(models => [$model]);
                $cleanup_succeeded->execute;
                $cleanup_rv{$cleanup_succeeded->result}++;
                $track_change->();
            }
        }
        elsif ($latest_build_status eq 'Succeeded') {
            $action = 'none';
        }
        elsif (should_review_model($model)) {
            $action = 'review';
        }
        else {
            $action = 'rebuild';
            if ($self->auto) {
                $request_build->();
                $track_change->();
            }
        }

        my $has_hidden_status = grep { lc $_ eq lc $latest_build_status } @hide_statuses;
        my $should_hide_none = $self->hide_no_action_needed && $action eq 'none';
        unless ($has_hidden_status || $should_hide_none) {
            $self->print_message(join "\t", $model_id, $action, $latest_build_status, $first_nondone_step, $latest_build_revision, $model_name, $pp_name, $fail_count);
        }

        # help avoid bad state in the larger auto_batch_size transaction
        unless ($per_model_txn->commit) {
            $per_model_txn->rollback;
        }

        if ($change_count > $self->auto_batch_size) {
            $commit->();
        }
    }

    $self->print_message("Requested builds for $build_requested_count models.") if $build_requested_count;
    $self->print_message("Cleaned up " . $cleanup_rv{1} . ".") if $cleanup_rv{1};
    $self->print_message("Failed to clean up " . $cleanup_rv{0} . ".") if $cleanup_rv{0};
    return 1;
}


sub print_message {
    my $self = shift;
    my $msg = shift;
    print STDOUT $msg . "\n";
    return 1;
}


sub should_review_model {
    my $model = shift;

    # If the latest build succeeded then we're happy.
    return if latest_build_succeeded($model);

    my @builds = $model->builds;
    return if @builds < 1;

    my $latest_status = $model->latest_build->status;
    return 1 if $latest_status eq 'Unstartable';

    return if @builds == 1;

    # If it has failed >3 times in a row then submit for review.
    return 1 if model_has_failed_to_many_times($model);

    # If it hasn't made progress since last time then submit for review.
    return 1 unless model_has_progressed($model);

    return;
}


sub latest_build_succeeded {
    my $model = shift;

    my $latest_build = $model->latest_build;
    return ($latest_build->status eq 'Succeeded');
}


my $max_fails = 5;
sub model_has_failed_to_many_times {
    my $model = shift;

    my @builds = reverse $model->builds;
    return unless @builds;

    my $n = (@builds >= $max_fails ? ($max_fails - 1) : $#builds);
    my @last_n_builds = @builds[0..$n];

    my @failed_builds = grep { $_->status eq 'Failed' } @last_n_builds;
    return (@failed_builds == $max_fails );
}


sub model_has_progressed {
    my $model = shift;

    my @builds = $model->builds;
    die unless (@builds > 1);

    my $latest_build = $model->latest_build;
    my %latest_status = workflow_status($latest_build);

    my $previous_build = previous_build($latest_build);
    my %previous_status = workflow_status($previous_build);

    my $status_are_different = status_compare(\%latest_status, \%previous_status);

    return $status_are_different;
}


sub latest_build_is_succeeded {
    my $model = shift;

    my $latest_build = $model->latest_build;

    return (
        $latest_build
        && $latest_build->status
        && $latest_build->status eq 'Succeeded'
    );
}


sub previous_build {
    my $build = shift;

    my $model = $build->model;
    my @prior_builds = grep { $_->date_scheduled lt $build->date_scheduled } $model->builds;

    my $previous_build = @prior_builds ? $prior_builds[-1] : undef;
    return $previous_build;
}


sub workflow_status {
    my $build = shift;

    # The following is wrapped in a transaction and try to protect against
    # "corrupt" Workflow::Models.
    my $tx = UR::Context::Transaction->begin();
    my %status = try {
        my $build_instance = $build->newest_workflow_instance;
        my @child_instances = $build_instance ? $build_instance->sorted_child_instances : ();

        my %status;
        for my $child_instance (@child_instances) {
            (my $name = $child_instance->name) =~ s/^[0-9]+\s+//;
            $status{$name} = $child_instance->status;
        }
        $tx->commit();
        return %status;
    } catch {
        $tx->rollback();
        return;
    };

    return %status;
}


sub status_compare { # http://stackoverflow.com/q/540229
    my %a = %{ shift @_ };
    my %b = %{ shift @_ };
    for my $key (keys %a) {
          return 1 unless defined $b{$key} and $a{$key} eq $b{$key};
          delete $a{$key};
          delete $b{$key};
    }
    return 1 if keys %b;
    return 0;
}


sub find_first_nondone_step {
    my $build = shift;

    # The following is wrapped in a transaction and try to protect against
    # "corrupt" Workflow::Models.
    my $tx = UR::Context::Transaction->begin();
    my $first_nondone_step = try {
        my $wf = $build->newest_workflow_instance;
        my $step =  _find_first_nondone_step_impl($wf);
        $tx->commit();
        return $step;
    } catch {
        $tx->rollback();
        return;
    };

    return $first_nondone_step;
}


sub _find_first_nondone_step_impl {
    my $parent_workflow_instance = shift;
    my @child_workflow_instances = $parent_workflow_instance->related_instances;

    my $failed_step;
    for my $child_workflow_instance (@child_workflow_instances) {
        $failed_step = _find_first_nondone_step_impl($child_workflow_instance);
        last if $failed_step;
    }
    # detect-variants is skipped because of the way the DV2 dispatcher works, not sure if this will work in general though
    my $name = $parent_workflow_instance->name;
    my $is_detect_variants = ($name =~ /^detect-variants/ || $name =~/^Detect\ Variants/);
    if ($parent_workflow_instance->status ne 'done' and not $failed_step and !$is_detect_variants) {
        $failed_step = $parent_workflow_instance->name;
    }

    return $failed_step;
}
