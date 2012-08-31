package Genome::Model::Tools::Smalt;

use strict;
use warnings;

use Genome;

my $DEFAULT = '0.5.5';

class Genome::Model::Tools::Smalt {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of smalt to use, default is $DEFAULT" },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run smalt or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools smalt ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}

my %VERSIONS = (
    '0.5.5'    => '/gsc/bin/smalt_x86_64',
    '0.2.8'    => '/gscmnt/sata820/info/medseq/alignment-test/smalt/smalt-0.2.8/smalt_x86-64',
    'smalt' => 'smalt',
);

sub smalt_path {
    my $self = $_[0];
    return $self->path_for_smalt_version($self->use_version);
}

sub available_smalt_versions {
    my $self = shift;
    return keys %VERSIONS;
}

sub path_for_smalt_version {
    my $class = shift;
    my $version = shift;

    if (defined $VERSIONS{$version}) {
        return $VERSIONS{$version};
    }
    die('No path for smalt version '. $version);
}

sub default_smalt_version {
    die "default samtools version: $DEFAULT is not valid" unless $VERSIONS{$DEFAULT};
    return $DEFAULT;
}

sub default_version { return default_smalt_version; }
        

1;

