package Genome::Model::Tools::Shrimp2;

use strict;
use warnings;

use Genome;

my $DEFAULT = '2.0.1';

class Genome::Model::Tools::Shrimp2 {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of shrimp2 to use, default is $DEFAULT" },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run SHRiMP2 or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools shrimp2 ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}

my %VERSIONS = (
#    '2.0.1'   => $ENV{GENOME_SW} . '/shrimp2/SHRiMP_2_0_1/bin/gmapper',
    '2.0.1'   => '/gscmnt/sata820/info/medseq/alignment-test/shrimp2/SHRiMP_2_0_1/bin/gmapper',
    'shrimp2' => 'shrimp2',
);

sub shrimp2_path {
    my $self = $_[0];
    return $self->path_for_shrimp2_version($self->use_version);
}

sub available_shrimp2_versions {
    my $self = shift;
    return keys %VERSIONS;
}

sub path_for_shrimp2_version {
    my $class = shift;
    my $version = shift;

    if (defined $VERSIONS{$version}) {
        return $VERSIONS{$version};
    }
    die('No path for shrimp2 version '. $version);
}

sub default_shrimp2_version {
    die "default shrimp2 version: $DEFAULT is not valid" unless $VERSIONS{$DEFAULT};
    return $DEFAULT;
}
        

1;

