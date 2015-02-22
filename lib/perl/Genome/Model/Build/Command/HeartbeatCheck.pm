package Genome::Model::Build::Command::HeartbeatCheck;

use strict;
use warnings;

class Genome::Model::Build::Command::HeartbeatCheck {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Builds to check heartbeat.',
        },
        verbose => {
            is => 'Boolean',
            default => 1,
            doc => 'Shows reason for heartbeat failure.',
        },
    ],
};

sub help_detail {
    "Check supplied builds' heartbeat."
}

sub execute {
    my $self = shift;

    my @builds = $self->builds;
    for my $build (@builds) {
        $build->dump_status_messages(1);
        my $heartbeat = ($build->heartbeat(verbose => $self->verbose) || '0');
        print join("\t", $build->id, $heartbeat) . "\n";
    }

    return 1;
}

1;

