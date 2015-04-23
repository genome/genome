package Genome::Model::Tools::SpeedSeq::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SpeedSeq::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has_param => [
        version => {
            is => 'Version',
            doc => 'SpeedSeq version to be used',
            is_optional => 1,
        },
        config_file => {
            is => 'Text',
            doc => 'path to speedseq.config file (default: same directory as speedseq)',
            is_optional => 1,
        },
        verbose => {
            is => 'Boolean',
            doc => 'verbose',
            is_optional => 1,
        },
    ],
};

sub help_detail {
    "A flexible framework for rapid genome analysis and interpretation"
}

my %SPEEDSEQ_VERSIONS = (
    'test'     => '/gscmnt/gc2719/halllab/bin/speedseq',
);

sub available_versions {
    my $self = shift;
    return keys(%SPEEDSEQ_VERSIONS);
}

sub path_for_version {
    my ($class, $version) = @_;
    unless ($version) {
        die('No version defined!');
    }
    my $path = $SPEEDSEQ_VERSIONS{$version};
    return $path if defined $path;
    die 'No path found for version: '. $version;
}

sub speedseq_path {
    my $self = shift;
    return $self->path_for_version($self->version);
}

1;
