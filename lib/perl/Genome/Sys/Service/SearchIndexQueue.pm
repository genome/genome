package Genome::Sys::Service::SearchIndexQueue;
use strict;
use warnings;
use Genome;

class Genome::Sys::Service::SearchIndexQueue {
    is => ['Genome::Sys::Service','UR::Singleton'],
    doc => "SearchIndexQueue service"
};

sub host { "solr.gsc.wustl.edu" }

sub restart_command { "service genome-search-index-queue restart" }

sub stop_command { "service genome-search-index-queue stop" }

sub log_path { "/var/log/genome/search-index-queue.log" }

sub pid_name { "genome sys search index daemon" }

sub url { }
