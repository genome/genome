package Genome::Sys::Service::Command::Restart;
use warnings;
use strict;
use Genome;

class Genome::Sys::Service::Command::Restart {
    is => 'Command::V2',
    has => [
        service => {
            is => 'Genome::Sys::Service',
            doc => 'The service which should be restarted',
            shell_args_position => 1,
        },
        output_basedir => {
            is => 'String',
            doc => 'Directory to put stdout/stderr',
            is_optional => 1,
        }
    ],
};

sub execute {

    my $self = shift;
    my $service = $self->service();

    my $cmd = join(' ', 'sudo ssh', $service->host(), $service->restart_command());
#    my $e = system($cmd);
#    die 'Error with cmd: $cmd' if $e;

    warn "would have executed cmd: $cmd";#
    return 1;
}


