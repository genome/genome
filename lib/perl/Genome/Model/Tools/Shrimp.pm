package Genome::Model::Tools::Shrimp;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Shrimp {
    is => 'Command',
    has_input => [
        use_version => {
            is => 'Version',
            default_value => '1.0.2',
            doc => "Version of shrimp to use(default_value=1.0.2)"
        },
        read_space => {
            is => 'Text',
            doc => "'ls', 'cs', or 'hs' for letter-space (454, Illumina/Solexa), colour-space (AB SOLiD), and Helicos-space (Helicos HeliScope SMS 2-pass reads), respectively.(default_value=ls)",
            default_value => 'ls',
            valid_values => ['ls','cs','hs'],
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run SHRiMP or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools shrimp ...    
EOS
}

sub help_detail {                           
    return <<EOS 
More information about the SHRiMP suite of tools can be found at http://compbio.cs.toronto.edu/shrimp/.
EOS
}

sub shrimp_path {
    my $self = $_[0];
    return $self->path_for_shrimp_version($self->use_version);
}

my %SHRIMP_VERSIONS = (
		    '1.0.2' => $ENV{GENOME_SW} . '/shrimp/SHRiMP_1_0_2/bin',
                );

sub available_shrimp_versions {
    my $self = shift;
    return keys %SHRIMP_VERSIONS;
}

sub path_for_shrimp_version {
    my $class = shift;
    my $version = shift;

    if (defined $SHRIMP_VERSIONS{$version}) {
        return $SHRIMP_VERSIONS{$version};
    }
    die('No path for shrimp version '. $version);
}


1;

