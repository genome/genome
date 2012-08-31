package Genome::Model::Tools::Blat;

use strict;
use warnings;

use Genome;
use File::Basename;


#declare a default version here
##########################################
my $DEFAULT = '34';

class Genome::Model::Tools::Blat {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of Blat to use, default is $DEFAULT" },
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
    "Tools to run Blat or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools blat ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}


my %BLAT_VERSIONS = (
	'34' => '/gsc/bin/blat',
    'blat'   => 'blat',
);


sub blat_path {
    my $self = $_[0];
    return $self->path_for_blat_version($self->use_version);
}

sub available_blat_versions {
    my $self = shift;
    return keys %BLAT_VERSIONS;
}

sub path_for_blat_version {
    my $class = shift;
    my $version = shift;

    if (defined $BLAT_VERSIONS{$version}) {
        return $BLAT_VERSIONS{$version};
    }
    die('No path for Blat version '. $version);
}

sub default_blat_version {
    die "default samtools version: $DEFAULT is not valid" unless $BLAT_VERSIONS{$DEFAULT};
    return $DEFAULT;
}
        

1;

