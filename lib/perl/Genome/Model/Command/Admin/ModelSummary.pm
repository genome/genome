package Genome::Model::Command::Admin::ModelSummary;

class Genome::Model::Command::Admin::ModelSummary {
    is => 'Genome::Command::Base',
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
        auto => {
            is => 'Boolean',
            default => 0,
            doc => 'auto act based on recommended actions',
        },
        auto_batch_size  => {
            is => 'Integer',
            default => 0,
            doc => 'When "--auto" is true, commit() to the database after this many Models have been identified for cleanup or restart.  0 means wait until the end to commit everything',
        }
    ],
};

use strict;
use warnings;
use Genome;

sub execute {
    my $self = shift;

    if ($self->auto_batch_size > 0 and !$self->auto) {
        $self->error_message("Cannot use --auto_batch_size unless --auto is true");
        return;
    }

    my @models = $self->models;
    my @hide_statuses = $self->hide_statuses;

    #preload data
    my @builds = Genome::Model::Build->get(
        model_id => [map($_->id, @models)],
        -hint => ['the_master_event'],
    );

    my $synchronous = ($self->auto and $self->auto_batch_size);  # whether we should start builds as we go or wait until the end

    my $models_changed_count = 0;
    my $build_requested_count = 0;

    my %cleanup_rv;
    my @models_to_cleanup;
    my $cleanup_build = sub {
        my $model = shift;
        if ($synchronous) {
            my $cmd = Genome::Model::Command::Admin::CleanupSucceeded->create(models => [$model]);
            $cleanup_rv{$cmd->execute}++;
            if ($self->auto_batch_size and ++$models_changed_count > $self->auto_batch_size) {
                $self->status_message("\n\nCommitting progress so far...");
                die "Commit failed!  Bailing out." unless UR::Context->commit;
                $self->status_message("Continuing...\n");
                $models_changed_count = 0;
            }
        } else {
            push @models_to_cleanup, $model;
        }
    };

    # Header for the report produced at the end of the loop
    $self->print_message(join("\t", qw(model_id action latest_build_status first_nondone_step latest_build_rev model_name pp_name fail_count)));
    $self->print_message(join("\t", qw(-------- ------ ------------------- ------------------ ---------------- ---------- ------- ----------)));

    EXAMIME_MODELS:
    for my $model (@models) {
        my $latest_build        = ($model        ? $model->latest_build  : undef);
        my $latest_build_status = ($latest_build ? $latest_build->status : '-');
        $latest_build_status = 'Requested' if $model->build_requested;

        my $fail_count   = ($model ? scalar $model->failed_builds     : undef);
        my $model_id     = ($model ? $model->id                       : '-');
        my $model_name   = ($model ? $model->name                     : '-');
        my $pp_name      = ($model ? $model->processing_profile->name : '-');

        my $first_nondone_step = '-';
        eval {
            if ($latest_build) {
                my $parent_workflow_instance = $latest_build->newest_workflow_instance;
                $first_nondone_step = find_first_nondone_step($parent_workflow_instance) || '-';
            }
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
        elsif ($latest_build->status eq 'Scheduled' || $latest_build->status eq 'Running' || $model->build_requested) {
            $action = 'none';
        }
        elsif ($latest_build && $latest_build->status eq 'Succeeded') {
            $action = 'cleanup';
            $cleanup_build->($model);
        }
        elsif (should_review_model($model)) {
            $action = 'review';
        }
        else {
            $action = 'rebuild';
            $model->build_requested(1, 'ModelSummary recommends rebuild');
            $build_requested_count++;
        }

        next if (grep { lc $_ eq lc $latest_build_status } @hide_statuses);
        $self->print_message(join "\t", $model_id, $action, $latest_build_status, $first_nondone_step, $latest_build_revision, $model_name, $pp_name, $fail_count);
    }

    $synchronous = 1;  # now we'll actually cleanup models if they were just being noted before
    if ($self->auto and @models_to_cleanup) {
        for my $model (@models_to_cleanup) {
            $cleanup_build->($model);
        }
    }

    $self->print_message("Requested builds for $build_requested_count models") if $build_requested_count;
    $self->print_message("Cleaned up " . $cleanup_rv{1}) if $cleanup_rv{1};
    $self->print_message("Failed to cleanup " . $cleanup_rv{0}) if $cleanup_rv{0};
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
    my @prior_builds = grep { $_->id < $build->id } $model->builds;

    my $previous_build = @prior_builds ? $prior_builds[-1] : undef;
    return $previous_build;
}


sub workflow_status {
    my $build = shift;

    my %status = eval {
        my $build_instance = $build->newest_workflow_instance;
        my @child_instances = $build_instance->sorted_child_instances if $build_instance;

        my %status;
        for my $child_instance (@child_instances) {
            (my $name = $child_instance->name) =~ s/^[0-9]+\s+//;
            $status{$name} = $child_instance->status;
        }
        return %status;
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
    my $parent_workflow_instance = shift;
    my @child_workflow_instances = $parent_workflow_instance->related_instances;

    my $failed_step;
    for my $child_workflow_instance (@child_workflow_instances) {
        $failed_step = find_first_nondone_step($child_workflow_instance);
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
