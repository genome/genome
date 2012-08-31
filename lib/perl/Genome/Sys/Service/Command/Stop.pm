package Genome::Sys::Service::Command::Stop;
use warnings;
use strict;
use Genome;

class Genome::Sys::Service::Command::Stop {
    is => 'Command::V2',
    has => [
        service => {
            is => 'Genome::Sys::Service',
            doc => 'The service which should be stopped',
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

    my $cmd = join(' ', 'ssh', $service->host(), $service->stop_command());
#    my $e = system($cmd);
#    die 'Error with cmd: $cmd' if $e;

    warn "would have executed cmd: $cmd";#

    return 1;
}


