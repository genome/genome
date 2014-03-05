package Genome::File::Fetch;

use strict;
use warnings;

use Carp qw(croak);
use Exporter qw(import);
use File::Fetch qw();
use Genome::Sys qw();

our @EXPORT_OK = qw(fetch);

sub fetch {
    my ($uri, $to) = @_;

    unless (defined $to) {
        $to = Genome::Sys->create_temp_directory;
    }

    local $File::Fetch::TIMEOUT = 60;

    my $ff = File::Fetch->new(uri => $uri);

    my $where = $ff->fetch(to => $to);
    unless($where) {
        my $err = $ff->error() || "unknown error while fetching $uri";
        croak $err;
    }

    return $where;
}

1;
