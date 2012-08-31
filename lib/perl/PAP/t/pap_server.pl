#!/usr/bin/env genome-perl

use strict;

use above 'Workflow::Server::Hub';
use Workflow::Server::HTTPD;
use Workflow::Server::UR;
use PAP;

POE::Kernel->stop();

my $pid = fork;
if ($pid) {
    print "$$ parent\n";
    Workflow::Server::Hub->start;
} elsif (defined $pid) {
    print "$$ child\n";

    Workflow::Server::HTTPD->start;
    Workflow::Server::UR->start;
} else {
    warn "no child?";
}

