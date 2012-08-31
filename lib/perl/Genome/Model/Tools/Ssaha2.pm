package Genome::Model::Tools::Ssaha2;

use strict;
use warnings;

use Genome;

my $DEFAULT = '2.5';

class Genome::Model::Tools::Ssaha2 {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of ssaha2 to use, default is $DEFAULT" },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run SSAHA2 or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools ssaha2 ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}

my %VERSIONS = (
    '2.5'    => $ENV{GENOME_SW} . '/ssaha2/ssaha2_v2.5_x86_64/ssaha2',
    'ssaha2' => 'ssaha2',
);

sub ssaha2_path {
    my $self = $_[0];
    return $self->path_for_ssaha2_version($self->use_version);
}

sub available_ssaha2_versions {
    my $self = shift;
    return keys %VERSIONS;
}

sub path_for_ssaha2_version {
    my $class = shift;
    my $version = shift;

    if (defined $VERSIONS{$version}) {
        return $VERSIONS{$version};
    }
    die('No path for ssaha2 version '. $version);
}

sub default_ssaha2_version {
    die "default samtools version: $DEFAULT is not valid" unless $VERSIONS{$DEFAULT};
    return $DEFAULT;
}

sub default_version { return default_ssaha2_version; }
        

1;

