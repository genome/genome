package Genome::Sys::Service::Web;
use strict;
use warnings;
use Genome;

class Genome::Sys::Service::Web {
    is => ['Genome::Sys::Service','UR::Singleton'],
    doc => "Web service"
};

sub host { "apipe.gsc.wustl.edu" }

sub restart_command { "/etc/init.d/genome_view restart" }

sub stop_command { "/etc/init.d/genome_view stop" }

sub log_path { "/var/log/kom/genome_view.log" }

sub pid_name { "perl-fcgi-pm" }

sub url { $ENV{GENOME_SYS_SERVICES_SEARCH_URL} }
