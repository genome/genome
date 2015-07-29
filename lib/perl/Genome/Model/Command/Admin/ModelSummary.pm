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
    has_optional => [
        max_fail_count => {
            is => 'Number',
            default => 5,
            doc => 'Rebuild will not be recommended for models that have failed this many times or more. (If undefined, any fail count will be allowed.)',
        },
    ],
};


sub execute {
    my $self = shift;

    my @models = $self->models;
    my @hide_statuses = $self->hide_statuses;

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

        my ($latest_build, $latest_build_status) = $self->_build_and_status_for_model($model);

        my $failure_set  = $self->failure_build_set($model);
        my $fail_count   = $failure_set->count;
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
        }
        elsif ($latest_build_status eq 'Scheduled' || $latest_build_status eq 'Running' || $latest_build_status eq 'Requested') {
            $action = 'none';
        }
        elsif ($latest_build && $latest_build_status eq 'Succeeded' && $fail_count) {
            $action = 'cleanup';
        }
        elsif ($latest_build_status eq 'Succeeded') {
            $action = 'none';
        }
        elsif ($self->should_review_model($model, $latest_build_status)) {
            $action = 'review';
        }
        else {
            $action = 'rebuild';
        }

        #take action if requested
        if ($self->auto) {
            if ($action eq 'rebuild' or $action eq 'build-needed') {
                $request_build->();
                $track_change->();
            }
            elsif ($action eq 'cleanup') {
                my $cleanup_succeeded = Genome::Model::Command::Admin::CleanupSucceeded->create(models => [$model]);
                $cleanup_succeeded->execute;
                $cleanup_rv{$cleanup_succeeded->result}++;
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

sub _build_and_status_for_model {
    my $self = shift;
    my $model = shift;

    my $build_iterator = $model->build_iterator(
        'status not like' => 'Abandoned',
        '-order_by' => '-created_at',
    );
    my $latest_build        = $build_iterator->next;
    my $latest_build_status = ($latest_build ? $latest_build->status : '-');
    $latest_build_status = 'Requested' if $model->build_requested;

    return ($latest_build, $latest_build_status);
}


sub print_message {
    my $self = shift;
    my $msg = shift;
    print STDOUT $msg . "\n";
    return 1;
}


sub should_review_model {
    my $self = shift;
    my $model = shift;
    my $latest_status = shift;

    # If the latest build succeeded then we're happy.
    return if $latest_status eq 'Succeeded';

    # If there are no builds, nothing to review.
    return if $latest_status eq '-';

    return 1 if $latest_status eq 'Unstartable';

    # If it has failed >X times in a row then submit for review.
    return 1 if $self->model_has_failed_too_many_times($model);

    # If it hasn't made progress since last time then submit for review.
    return 1 unless $self->model_has_progressed($model);

    return;
}


sub model_has_failed_too_many_times {
    my $self = shift;
    my $model = shift;

    my $max_fails = $self->max_fail_count;
    return unless defined $max_fails;

    my $failure_set = $self->failure_build_set($model);

    return ($failure_set->count >= $max_fails);
}


sub model_has_progressed {
    my $self = shift;
    my $model = shift;

    my $failure_set = $self->failure_build_set($model);

    #a first build has always made progress
    return 1 unless ($failure_set->count > 1);

    my $it = $failure_set->member_iterator(-order_by => ['-created_at']);
    my $latest_build = $it->next;
    my $latest_error = determine_error_for_build($latest_build);
    return unless $latest_error;

    my $previous_build = $it->next;
    my $previous_error = determine_error_for_build($previous_build);
    return unless $previous_error;

    return $latest_error ne $previous_error;
}

sub failure_build_set {
    my $self = shift;
    my $model = shift;

    return Genome::Model::Build->define_set(
        model_id => $model->id,
        status => ['Failed', 'Unstartable', 'Unknown'],
    );
}


sub determine_error_for_build {
    my $build = shift;

    return unless ($build->status eq 'Failed' or $build->status eq 'Unstartable');

    my $cmd = Genome::Model::Build::Command::DetermineError->execute(
        build => $build,
        display_results => 0,
    );

    return $cmd->get_failed_key;
}


sub find_first_nondone_step {
    my $build = shift;

    return if grep { $_ eq $build->status } ('Succeeded', 'Scheduled', 'Unstartable', 'New');

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

1;
