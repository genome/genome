package Genome::Sys::Service::Memcache;
use strict;
use warnings;
use Genome;

class Genome::Sys::Service::Memcache {
    is => ['Genome::Sys::Service','UR::Singleton'],
    doc => "Memcache service"
};

{
    my $host;
    sub host {
        unless (defined $host) {
            $host = $ENV{GENOME_SYS_SERVICES_MEMCACHE};
            $host =~ s/:\w+$//;  # Strip out port designation
        }
        return $host;
    }
}

sub restart_command { "/etc/init.d/memcached restart" }

sub stop_command { "/etc/init.d/memcached stop" }

sub log_path { "/var/log/memcached.log" }

sub pid_name { "memcached" }

sub url {}

sub status {
    if(Genome::Memcache->server->stats->{total}) {
        return 'running';
    } else {
        return 'stopped';
    }
}
