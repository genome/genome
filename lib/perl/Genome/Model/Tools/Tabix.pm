package Genome::Model::Tools::Tabix;

use strict;
use warnings;

use Genome;
use Carp qw/confess/;
use Sys::Hostname;

class Genome::Model::Tools::Tabix {
    is  => 'Command',
    is_abstract => 1,
};


sub help_brief {
    "Tools to run tabix";
}

sub help_synopsis {
    "gmt tabix ...";
}

sub help_detail {                           
    "used to invoke tabix commands";
}

sub tabix_path {
    my ($self, $ver) = @_;
    my $path = "/usr/bin/tabix";

    my $hostname = Sys::Hostname::hostname();
    if (! -x $path) {
        confess "Failed to find executable tabix at $hostname:$path!";
    }
    return $path;
}

1;
