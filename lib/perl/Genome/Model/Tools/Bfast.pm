package Genome::Model::Tools::Bfast;

use strict;
use warnings;

use Genome;
use File::Basename;


#declare a default version here
##########################################
my $DEFAULT = '0.6.4d';

class Genome::Model::Tools::Bfast {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of Bfast to use, default is $DEFAULT" },
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
    "Tools to run Bfast or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools bfast ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}


my %BFAST_VERSIONS = (
	'0.6.4d' => '/gscmnt/sata820/info/medseq/alignment-test/bfast/bfast',
    'bfast'   => 'Bfast',
);


sub bfast_path {
    my $self = $_[0];
    return $self->path_for_bfast_version($self->use_version);
}

sub available_bfast_versions {
    my $self = shift;
    return keys %BFAST_VERSIONS;
}

sub path_for_bfast_version {
    my $class = shift;
    my $version = shift;

    if (defined $BFAST_VERSIONS{$version}) {
        return $BFAST_VERSIONS{$version};
    }
    die('No path for Bfast version '. $version);
}

sub default_bfast_version {
    die "default samtools version: $DEFAULT is not valid" unless $BFAST_VERSIONS{$DEFAULT};
    return $DEFAULT;
}
        

1;

