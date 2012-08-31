package Genome::Model::Tools::Crossmatch;
use strict;
use warnings;
use Genome;
use File::Basename;

my $DEFAULT = '1.080721';

class Genome::Model::Tools::Crossmatch {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of Crossmatch to use, default is $DEFAULT" },
        arch_os => {
                    calculate => q|
                            my $arch_os = `uname -m`;
                            chomp($arch_os);
                            return $arch_os;
                        |
                },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run Crossmatch or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools crossmatch ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}

my %CROSSMATCH_VERSIONS = (
	'1.080721' => '/gsc/bin/cross_match',
	'test' => '/gsc/bin/cross_match.test',
    'crossmatch'   => 'cross_match',
);

sub crossmatch_path {
    my $self = $_[0];
    return $self->path_for_crossmatch_version($self->use_version);
}

sub available_crossmatch_versions {
    my $self = shift;
    return keys %CROSSMATCH_VERSIONS;
}

sub path_for_crossmatch_version {
    my $class = shift;
    my $version = shift;

    if (defined $CROSSMATCH_VERSIONS{$version}) {
        return $CROSSMATCH_VERSIONS{$version};
    }
    die('No path for Crossmatch version '. $version);
}

sub default_crossmatch_version {
    die "default samtools version: $DEFAULT is not valid" unless $CROSSMATCH_VERSIONS{$DEFAULT};
    return $DEFAULT;
}

sub default_version { return default_crossmatch_version; }

1;

