package Genome::Model::Command::Admin::CleanupSucceeded;

class Genome::Model::Command::Admin::CleanupSucceeded {
    is => 'Genome::Command::Base',
    doc => 'Abandon failed builds for models with latest build succeeded.',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            is_optional => 1,
            shell_args_position => 1,
            doc => 'Models to check; resolved by Genome::Command::Base',
        },
    ],
};

use strict;
use warnings;
use Genome;

sub help_detail {
    'Abandon failed builds for models with latest build succeeded.'
}

sub execute {
    my $self = shift;

    my $user = Genome::Sys->username;
    my @builds_to_abandon;
    for my $model ($self->models) {
        my $latest_build = $model->latest_build;
        next unless $latest_build;
        next unless ($latest_build->status eq 'Succeeded');
        my @builds = $model->builds;
        for my $build (@builds) {
            next if ($build->id eq $latest_build->id);
            my $status = $build->status;
            next unless ($status eq 'Failed' or $status eq 'Unstartable' or $status eq 'New');
            next if ($build->run_by ne $user);
            push @builds_to_abandon, $build;
        }
    }

    if (@builds_to_abandon) {
        my $abandon_failed = Genome::Model::Build::Command::Abandon->create(
            builds => \@builds_to_abandon,
        );
        $abandon_failed->dump_status_messages(1);
        return $abandon_failed->execute();
    }

    return 1;
}
