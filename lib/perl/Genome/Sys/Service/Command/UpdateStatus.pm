package Genome::Sys::Service::Command::UpdateStatus;
use warnings;
use strict;
use Genome;

use LWP::UserAgent;

class Genome::Sys::Service::Command::UpdateStatus{
    is => 'Command::V2',
    has => [
        service => {
            is => 'Genome::Sys::Service',
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

    my ($self) = @_;

    $self->update_status();
    $self->update_pid_status();

    return 1;
}

sub update_status {

    my ($self) = @_;

    my $service = $self->service();

    return if !$service->url();

    my $ua = LWP::UserAgent->new();
    $ua->cookie_jar({ file => "/tmp/nomnomcookies.txt" });
    $ua->timeout(10);
    my $r = $ua->get($service->url());

    my $old_status = $service->status();
    my $status;

    if ($r->is_success) {
        $status = 'running';
    } else {
        $status = 'stopped';
    }

    if ($status ne $old_status) {
        $service->status($status);
        print 'Service "'. $service->name() . "\" status: $status (was $old_status)\n";
    }

    return 1;
}

sub update_pid_status {

    my ($self) = @_;

    my $service = $self->service();
    my $pid_name = $service->pid_name();

    my $cmd = 'ssh ' . $service->host() . " \"ps aux | grep $pid_name | grep -v grep\"";
    warn $cmd;
    my $o = `$cmd`;

    my $old_pid_status = $service->pid_status();
    my $pid_status;
    if ($o) {
        $pid_status = 'running';
    } else {
        $pid_status = 'stopped';
    }

    if ($pid_status ne $old_pid_status) {
        $service->pid_status($pid_status);
        print 'Service "'. $service->name() . " status: $pid_status (was $old_pid_status)\n";
    }
 
    return 1;
}


