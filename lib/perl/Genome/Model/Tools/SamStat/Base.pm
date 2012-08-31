package Genome::Model::Tools::SamStat::Base;

use strict;
use warnings;

use Genome;

my $SAMSTAT_DEFAULT = '1.08';

class Genome::Model::Tools::SamStat::Base {
	is => 'Command::V2',
	is_abstract => 1,
        has_input => [
            use_version => {
                is  => 'Version',
                doc => 'samstat version to be used.',
                is_optional   => 1,
                default_value => $SAMSTAT_DEFAULT,
            },
        ],
};

my @SAMSTAT_VERSIONS = (
     '1.08' => $ENV{GENOME_SW} . '/samstat/samstat-1.08/src/samstat',
     '1.03' => $ENV{GENOME_SW} . '/samstat/samstat-1.03/src/samstat',
);

my %SAMSTAT_VERSIONS = @SAMSTAT_VERSIONS;

sub latest_version { $SAMSTAT_VERSIONS[0] }

sub path_for_samstat_version {
    my ($class, $version) = @_;
    $version ||= $SAMSTAT_DEFAULT;
    my $path = $SAMSTAT_VERSIONS{$version};
    return $path if defined $path;
    die 'No path found for samstat version: '.$version;
}

sub default_samstat_version {
    die "default samstat version: $SAMSTAT_DEFAULT is not valid" unless $SAMSTAT_VERSIONS{$SAMSTAT_DEFAULT};
    return $SAMSTAT_DEFAULT;
}

sub samstat_path {
    my $self = shift;
    return $self->path_for_samstat_version($self->use_version);
}

1;
