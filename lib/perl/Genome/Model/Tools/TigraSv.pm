package Genome::Model::Tools::TigraSv;

use strict;
use warnings;

use Genome; 
use File::Basename;
use POSIX;
use DateTime;
use IO::File;

my $DEFAULT = '20110321';

class Genome::Model::Tools::TigraSv {
    is  => 'Command',
    has => [
        use_version => { 
            is  => 'Version', 
            doc => "tigra_sv version to be used, default is $DEFAULT. ", 
            is_optional   => 1, 
            default_value => $DEFAULT,   
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run tigra_sv or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt tigra-sv ...    
EOS
}

sub help_detail {                           
    return <<EOS 
EOS
}


my %TIGRASV_VERSIONS = (
    20110321 => $ENV{GENOME_SW} . '/tigrasv/TIGRA_SV-20110321/tigra_sv',     
    '0.1'    => '/usr/bin/tigra-sv0.1',
);

sub available_tigrasv_versions {
    return keys(%TIGRASV_VERSIONS);
}

sub path_for_tigrasv_version {
    my ($class, $version) = @_;
    $version ||= $DEFAULT;
    my $path = $TIGRASV_VERSIONS{$version};
    return $path if defined $path;
    die 'No path found for tigra_sv version: '.$version;
}

sub default_tigrasv_version {
    die "default tigra_sv version: $DEFAULT is not valid" unless $TIGRASV_VERSIONS{$DEFAULT};
    return $DEFAULT;
}    
    
sub tigrasv_path {
    my $self = shift;
    return $self->path_for_tigrasv_version($self->use_version);
}



1;
