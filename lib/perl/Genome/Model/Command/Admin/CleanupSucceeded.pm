package Genome::Model::Command::Admin::CleanupSucceeded;

class Genome::Model::Command::Admin::CleanupSucceeded {
    is => 'Command::V2',
    roles => ['Genome::Model::Command::Submittable'],
    doc => 'Abandon unsuccessful builds that have been superceded.',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            is_optional => 1,
            shell_args_position => 1,
            doc => 'Models to check',
        },
        submit_abandon_jobs => {
            is => 'Boolean',
            default => 0,
            doc => 'if set, will use `genome analysis-project possible-volumes` to submit jobs for abandoning',
        },
    ],
};

use strict;
use warnings;
use Genome;

sub help_detail {
    'Abandon unsuccessful builds that have been superceded.'
}

sub execute {
    my $self = shift;

    my $user = Genome::Sys->username;
    my @builds_to_abandon;
    for my $m ($self->models) {
        my $latest_unabandoned_build = $m->builds(-order_by => ['-created_at'], -limit => 1, 'status !=' => 'Abandoned');
        next unless $latest_unabandoned_build;

        my %params = (
            'id !=' => $latest_unabandoned_build->id,
        );

        if($user ne 'apipe-builder') {
            $params{run_by} = $user;
        }

        # if we made it out of Unstartable then the Unstartable problem is
        # fixed so we can abandon previous Unstartable builds
        if ($latest_unabandoned_build->status ne 'Unstartable') {
            push @builds_to_abandon,
                $m->builds(%params, status => 'Unstartable');
        }

        # if we succeeded then previous failures can be abandoned
        if ($latest_unabandoned_build->status eq 'Succeeded') {
            push @builds_to_abandon,
                $m->builds(%params, 'status in' => [qw(Failed New)]);
        }
    }

    if (@builds_to_abandon) {
        $self->status_message(sprintf('Will abandon %d builds...',
            scalar(@builds_to_abandon)));

        if ($self->submit_abandon_jobs) {
            $self->_submit_abandon_jobs(@builds_to_abandon);
        } else {
            my $abandon_failed = Genome::Model::Build::Command::Abandon->create(
                builds => \@builds_to_abandon,
                show_display_command_summary_report => 0,
            );
            $abandon_failed->dump_status_messages(1);
            return $abandon_failed->execute();
        }
    }

    return 1;
}

sub _submit_abandon_jobs {
    my $self = shift;
    my @builds = @_;

    my %builds_by_anp;
    for my $build (@builds) {
        my $anp = $build->model->analysis_project;
        $builds_by_anp{$anp->id} //= [];
        push @{ $builds_by_anp{$anp->id} }, $build;
    }

    for my $anp_id (keys %builds_by_anp) {
        my $anp = Genome::Config::AnalysisProject->get($anp_id);

        $self->_submit_jobs (
            $anp,
            [qw(genome model build abandon), map $_->id, @{$builds_by_anp{$anp_id}}],
        );
    }

    return 1;
}

1;
