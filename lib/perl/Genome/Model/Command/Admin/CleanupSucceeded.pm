package Genome::Model::Command::Admin::CleanupSucceeded;

class Genome::Model::Command::Admin::CleanupSucceeded {
    is => 'Command::V2',
    doc => 'Abandon unsuccessful builds that have been superceded.',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            is_optional => 1,
            shell_args_position => 1,
            doc => 'Models to check',
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
        my $abandon_failed = Genome::Model::Build::Command::Abandon->create(
            builds => \@builds_to_abandon,
            show_display_command_summary_report => 0,
        );
        $abandon_failed->dump_status_messages(1);
        return $abandon_failed->execute();
    }

    return 1;
}
