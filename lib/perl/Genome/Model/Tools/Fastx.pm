package Genome::Model::Tools::Fastx;

use strict;
use warnings;

use Genome;
use File::Basename;

my $DEFAULT = '0.0.13';

class Genome::Model::Tools::Fastx {
    is => 'Command',
    has => [
        use_version => {
            is => 'Version',
            is_optional => 1,
            default_value => $DEFAULT,
            doc => "Version of fastx to use, default is $DEFAULT"
        },
    ],
    has_abstract_constant => [
        fastx_tool => { is => 'Text', },
    ],
};

sub help_brief {
    "Tools to run Fastx or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools fastx ...    
EOS
}

sub help_detail {
    return <<EOS
More information about the Fastx toolkit can be found at http://hannonlab.cshl.edu/fastx_toolkit/
EOS
}


my %FASTX_VERSIONS = (
    '0.0.13' => $ENV{GENOME_SW} . '/fastx/fastx_toolkit-0.0.13/bin',
    '0.0.10' => $ENV{GENOME_SW} . '/fastx/fastx_toolkit-0.0.10/bin',
    '0.0.7' => $ENV{GENOME_SW} . '/fastx/fastx_toolkit-0.0.7/bin',
);

sub fastx_tool_path {
    my $self = shift;
    unless ($self->fastx_tool) {
        return;
    }
    return $self->fastx_base_path .'/'. $self->fastx_tool;
}

sub fastx_base_path {
    my $self = $_[0];
    return $self->base_path_for_fastx_version($self->use_version);
}

sub available_fastx_versions {
    my $self = shift;
    return keys %FASTX_VERSIONS;
}

sub base_path_for_fastx_version {
    my $class = shift;
    my $version = shift;

    if (defined $FASTX_VERSIONS{$version}) {
        return $FASTX_VERSIONS{$version};
    }
    die('No path for fastx version '. $version);
}

sub default_fastx_version {
    die "default fastx version: $DEFAULT is not valid" unless $FASTX_VERSIONS{$DEFAULT};
    return $DEFAULT;
}

1;

