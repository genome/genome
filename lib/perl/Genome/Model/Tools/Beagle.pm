package Genome::Model::Tools::Beagle;

use strict;
use warnings;

use Genome;

my $DEFAULT_VERSION = "Test";

class Genome::Model::Tools::Beagle {
    is => ['Command'],
    doc => 'Tools and scripts to phase and impute data using Beagle software.',
    has => [
    version => {
        is => 'Text',
        default => $DEFAULT_VERSION,
    },
    ],
};

my %BEAGLE_VERSIONS = (
    'Test' => '/gscmnt/ams1161/info/model_data/kmeltzst/Software/b4.r1128.jar',
);

sub get_default_beagle_version {
    return $DEFAULT_VERSION;
}

sub path_for_version {
    my $class = shift;
    my $version = shift || $DEFAULT_VERSION;

    if($version eq 'latest') {
        return $class->path_for_latest_version;
    }

    unless(exists $BEAGLE_VERSIONS{$version}) {
        $class->error_message('No path found for Beagle Version ' . $version);
        die $class->error_message;
    }

    return $BEAGLE_VERSIONS{$version};
}
sub path_for_latest_version {
    my $class = shift;
    my $path = '/gscmnt/ams1161/info/model_data/kmeltzst/Software/b4.r1128.jar';

    unless(-e $path) {
        $class->error_message("Path to latest version not found $path!");
    }

    return $path;
}

sub default_version {
    my $class = shift;

    unless(exists $BEAGLE_VERSIONS{$DEFAULT_VERSION}) {
        $class->error_message('Default Beagle version (' . $DEFAULT_VERSION . ') is invalid.');
        die $class->error_message;
    }

    return $DEFAULT_VERSION;
}

sub available_varscan_versions {
    return keys(%BEAGLE_VERSIONS);
}

sub help_brief {
    "Tools and scripts to phase and impute data using Beagle software."
}

1;
