package Genome::Model::Tools::Novocraft;

use strict;
use warnings;

use Genome;
use File::Basename;

my $DEFAULT_VERSION = '2.05.20';

class Genome::Model::Tools::Novocraft {
    is => 'Command',
    has_param => [
        use_version => {
            is => 'Version',
            default_value => $DEFAULT_VERSION,
            doc => 'Version of novocraft to use. default_value='. $DEFAULT_VERSION,
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run novocraft or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools novocraft ...    
EOS
}

sub help_detail {                           
    return <<EOS 
More information about the novocraft suite of tools can be found at http://novocraft.sourceforege.net.
EOS
}

sub novoindex_path {
    my $self = shift;
    return $self->novocraft_path .'/novoindex';
}

sub novoalign_path {
    my $self = shift;
    return $self->novocraft_path .'/novoalign';
}

sub novocraft_path {
    my $self = $_[0];
    return $self->path_for_novocraft_version($self->use_version);
}
my %NOVOCRAFT_VERSIONS = (
                    '2.03.12' => $ENV{GENOME_SW} . '/novocraft/novocraft-2.03.12',
                    '2.04.02' => $ENV{GENOME_SW} . '/novocraft/novocraft-2.04.02',
                    '2.05.13' => $ENV{GENOME_SW} . '/novocraft/novocraft-2.05.13',
                    '2.05.20' => $ENV{GENOME_SW} . '/novocraft/novocraft-2.05.20',
                    '2.05.32' => $ENV{GENOME_SW} . '/novocraft/novocraft-2.05.33',
                    '2.07.11' => $ENV{GENOME_SW} . '/novocraft/novocraft-2.07.11', 
                    'novocraft'   => 'novoalign',
                );

sub available_novocraft_versions {
    my $self = shift;
    return keys %NOVOCRAFT_VERSIONS;
}

sub path_for_novocraft_version {
    my $class = shift;
    my $version = shift;

    if (defined $NOVOCRAFT_VERSIONS{$version}) {
        return $NOVOCRAFT_VERSIONS{$version} . '/novoalign';
    }
    die('No path for novocraft version '. $version);
}

sub path_for_novosam_version {
    my ($class, $version) = @_;
    $version ||= $DEFAULT_VERSION;
    my $path = $NOVOCRAFT_VERSIONS{$version} . '/novo2sam.pl';
    return $path if -s $path and -x $path;
    die "novosam path of version: $version is invalid\n";
}


sub default_novocraft_version {
    die "default novocraft version: $DEFAULT_VERSION is not valid" unless $NOVOCRAFT_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}

sub default_version { return default_novocraft_version; }

1;

