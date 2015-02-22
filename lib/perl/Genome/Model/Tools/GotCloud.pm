package Genome::Model::Tools::GotCloud;

use strict;
use warnings;

use Genome;
use Sys::Hostname;

class Genome::Model::Tools::GotCloud {
    is  => 'Command',
    is_abstract => 1,
};

sub help_brief {
    "Tools to run GotCloud.";
}

sub help_synopsis {
    "gmt gotcloud ...";
}

sub help_detail {
    "used to invoke gotcloud commands";
}

sub gotcloud_path {
    my ($self) = @_;
    my $path = "/usr/local/gotcloud/bin/gotcloud";
    my $hostname = Sys::Hostname::hostname();
    if (! -x $path) {
        die $self->error_message("Failed to find executable gotcloud version at $hostname:$path!");
    }
    return $path;
}

1;

