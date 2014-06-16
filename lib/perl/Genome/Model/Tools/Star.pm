package Genome::Model::Tools::Star;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use File::Temp;
use POSIX qw(floor);

my $DEFAULT_VERSION = '2.3.1z1';

class Genome::Model::Tools::Star {
    is => ['Command'],
    has_input => [
        use_version => {
            is      => 'string',
            doc     => 'version of star to use',
            default => $DEFAULT_VERSION,
        },
    ],
};


my %STAR_VERSIONS = (
    '2.3.1z1' => '/usr/bin/star2.3.1z1/STAR',
);

sub help_brief {
    "tools to work with star aligner"
}

sub help_detail {
    return <<EOS
EOS
}

sub star_path {
    my $self = shift;
    return $self->path_for_star_version($self->use_version);
}

sub available_star_versions {
    return keys %STAR_VERSIONS;
}

sub path_for_star_version {
    my ($self, $version) = @_;

    if (defined $STAR_VERSIONS{$version}) {
        my $path = $STAR_VERSIONS{$version};
        unless (-x $path) {
            die $self->error_message("star: $path is not executable");
        }
        return $path;
    }
    die 'No path for star version '. $version;
}

sub default_star_version {
    die "default star version: $DEFAULT_VERSION is not valid" unless $STAR_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

1;

