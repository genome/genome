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

    my @headers = qw(model_id action latest_build_status latest_build_rev model_name pp_name fail_count);

    # Header for the report produced at the end of the loop
    $self->print_message(join("\t", @headers));
    $self->print_message(join("\t", qw(-------- ------ ------------------- ---------------- ---------- ------- ----------)));

    my $build_requested_count = 0;
    my %cleanup_rv;
    my $abandon_count = 0;
    my $abandon_attempt_count = 0;
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

            $change_count = 0;
            $auto_batch_size_txn = UR::Context::Transaction->begin;
        } else {
            $auto_batch_size_txn->rollback;
        }
    };

    Genome::Model::Build::Set->class; #generate type before creating pool

    for my $model (@models) {
        my $guard = UR::Context::AutoUnloadPool->create();
        my $per_model_txn = UR::Context::Transaction->begin();

        my $summary = $self->generate_model_summary($model);

        #take action if requested
        my $action = $summary->{action};
        if ($self->auto) {
            if ($action eq 'rebuild' or $action eq 'build-needed') {
                $model->build_requested(1, 'ModelSummary recommended rebuild.');
                $build_requested_count++;
                $change_count++;
            }
            elsif ($action eq 'cleanup') {
                my $cleanup_succeeded = Genome::Model::Command::Admin::CleanupSucceeded->create(models => [$model]);
                $cleanup_succeeded->execute;
                $cleanup_rv{$cleanup_succeeded->result}++;
                $change_count++;
            } elsif ($action eq 'abandon') {
                my @builds = $self->failure_build_set($model)->members;
                local $ENV{UR_NO_REQUIRE_USER_VERIFY} = 1;
                my $abandon = Genome::Model::Build::Command::Abandon->create(
                    builds => \@builds,
                    header_text => 'Automatically abandoning builds for disabled model',
                    body_text => 'Triggered by ' . __PACKAGE__,
                    show_display_command_summary_report => 0,
                );
                $abandon->execute;
                $abandon_count += scalar(@builds) if $abandon->result;
                $abandon_attempt_count += scalar(@builds);
                $change_count++;
            }
        }

        my $has_hidden_status = grep { lc $_ eq lc $summary->{latest_build_status} } @hide_statuses;
        my $should_hide_none = $self->hide_no_action_needed && $action eq 'none';
        unless ($has_hidden_status || $should_hide_none) {
            $self->print_message(join "\t", @$summary{@headers});
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
    $self->print_message("Abandoned builds for ". $abandon_count . " (out of " . $abandon_attempt_count . " attempted).") if $abandon_attempt_count;
    return 1;
}

sub generate_model_summary {
    my $self = shift;
    my $model = shift;

    my %summary;

    my $failure_set  = $self->failure_build_set($model);
    $summary{fail_count}   = $failure_set->count;
    $summary{model_id}     = ($model ? $model->id                       : '-');
    $summary{model_name}   = ($model ? $model->name                     : '-');
    $summary{pp_name}      = ($model ? $model->processing_profile->name : '-');

    my ($latest_build, $latest_build_status) = $self->_build_and_status_for_model($model);
    $summary{latest_build_status} = $latest_build_status;

    my $latest_build_revision = $latest_build->software_revision if $latest_build;
    $latest_build_revision ||= '-';
    ($latest_build_revision) =~ /\/(genome-[^\/])/ if $latest_build_revision =~ /\/genome-[^\/]/;
    $summary{latest_build_rev} = $latest_build_revision;

    my $action;
    if ($latest_build_status eq 'Disabled') {
        $action = 'abandon';
    }
    elsif (!$latest_build && $model->build_requested){
        $action = 'none';
    }
    elsif (!$latest_build) {
        if($model->status eq 'Buildless') {
            $action = 'none';
        } else {
            $action = 'build-needed';
        }
    }
    elsif ($latest_build_status eq 'Scheduled' || $latest_build_status eq 'Running' || $latest_build_status eq 'Requested') {
        $action = 'none';
    }
    elsif ($latest_build && $latest_build_status eq 'Succeeded' && $summary{fail_count}) {
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
    $summary{action} = $action;

    return \%summary;
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
    $latest_build_status = 'Disabled' if $model->is_disabled;

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

    return 1 if $self->latest_failure_requires_review($model);

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


sub latest_failure_requires_review {
    my $self = shift;
    my $model = shift;

    my $failure_set = $self->failure_build_set($model);
    my $count = $failure_set->count;

    return unless $count > 0; #nothing to review

    my $it = $failure_set->member_iterator(-order_by => ['-created_at']);
    my $latest_build = $it->next;
    my $latest_error = determine_error_for_build($latest_build);
    return 1 unless $latest_error;

    return 1 if $self->_error_requires_review($latest_error);

    #if this error doesn't require review, retry first failures
    return unless $count > 1;

    my $previous_build = $it->next;
    my $previous_error = determine_error_for_build($previous_build);
    return 1 unless $previous_error;

    #review if consistently failing
    return $latest_error eq $previous_error;
}

sub _error_requires_review {
    my $self = shift;
    my $latest_error = shift;

    return 1 if ($latest_error =~ /TERM_MEMLIMIT/);
    return 1 if ($latest_error =~ /the allocation has been orphaned/);

    return;
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

1;
