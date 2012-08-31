package Genome::Model::Tools::Far::Base;

use strict;
use warnings;

use Genome;

my $FAR_DEFAULT     = '1.84';

class Genome::Model::Tools::Far::Base {
    is          => 'Command::V2',
    is_abstract => 1,
    has         => [
            use_version => {
                is  => 'Version',
                doc => "far version to be used, default is $FAR_DEFAULT. ",
                is_optional   => 1,
                default_value => $FAR_DEFAULT,
            },
    ],
};

my %FAR_VERSIONS = (
    '1.7'     => $ENV{GENOME_SW} . '/flexibleadapter/flexibleadapter-1.7/far',
    '1.84'    => $ENV{GENOME_SW} . '/flexibleadapter/flexibleadapter-1.84/build/far',
    '2.0'    => $ENV{GENOME_SW} . '/flexibleadapter/flexibleadapter-2.0/build/far',
    '2.17'    => '/usr/bin/far2.17.0',
);

sub available_far_versions {
    my $self = shift;
    return keys(%FAR_VERSIONS);
}

sub path_for_far_version {
    my ($class, $version) = @_;
    $version ||= $FAR_DEFAULT;
    my $path = $FAR_VERSIONS{$version};
    return $path if defined $path;
    die 'No path found for far version: '.$version;
}

sub default_far_version {
    die "default far version: $FAR_DEFAULT is not valid" unless $FAR_VERSIONS{$FAR_DEFAULT};
    return $FAR_DEFAULT;
}

sub far_path {
    my $self = shift;
    return $self->path_for_far_version($self->use_version);
}


1;
