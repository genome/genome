package Genome::Model::Tools::Squaredancer;

use warnings;
use strict;

use Genome;

my $DEFAULT_VERSION = '0.1';

class Genome::Model::Tools::Squaredancer{
    is => 'Command',
    has => [
        use_version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of squaredancer to use"
        },
     ],
};

my %SQUAREDANCER_VERSIONS = (
    '0.1' => '/gsc/scripts/opt/genome-stable/lib/perl/Genome/Model/Tools/Sv/SquareDancer.pl',
);


sub help_brief {
    "discovers structural variation using squaredancer",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
EOS
}

sub help_detail {                           
    return <<EOS 
This tool discovers structural variation for soft-clipped genomic or cDNA reads.
EOS
}

sub squaredancer_path {
    my $self = shift;
    return $self->path_for_squaredancer_version($self->use_version);
}

sub available_squaredancer_versions {
    my $self = shift;
    return keys %SQUAREDANCER_VERSIONS;
}

sub path_for_squaredancer_version {
    my ($self, $version) = @_;
    my $sd_path = $SQUAREDANCER_VERSIONS{$version};

    if (defined $sd_path) {
        unless (-s $sd_path and -x $sd_path) {
            $self->error_message("squaredancer executable path $sd_path for version $version is not valid");
            die $self->error_message;
        }
        return $sd_path;
    }
    die 'No path for squaredancer version '. $version;
}

sub default_squaredancer_version {
    die "default squaredancer version: $DEFAULT_VERSION is not valid" unless $SQUAREDANCER_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}
 
1;
