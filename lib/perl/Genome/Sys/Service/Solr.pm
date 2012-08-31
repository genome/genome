package Genome::Sys::Service::Solr;
use strict;
use warnings;
use Genome;

class Genome::Sys::Service::Solr {
    is => ['Genome::Sys::Service','UR::Singleton'],
    doc => "Solr service"
};

sub host { "solr.gsc.wustl.edu" }

sub restart_command { "/etc/init.d/tomcat6 restart" }

sub stop_command { "/etc/init.d/tomcat6 stop" }

sub log_path { "/var/log/tomcat6/catalina.log" }

sub pid_name { "tomcat6" }

sub url { "http://solr:8080/solr/admin" }
