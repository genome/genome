package Genome::Db::Tgi::Command::Util;

use strict;
use warnings;

use Exporter qw(import);
use File::Fetch qw();

our @EXPORT_OK = qw(fetch);

sub fetch {
    my ($uri, $to) = @_;
    local $File::Fetch::TIMEOUT = 60;
    my $ff = File::Fetch->new(uri => $uri);
    my $where = $ff->fetch(to => $to);
    unless($where) {
        my $err = $ff->error() || "unknown error while fetching $uri";
        die $err;
    }
    return $where;
}

1;
