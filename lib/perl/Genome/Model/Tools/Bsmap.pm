package Genome::Model::Tools::Bsmap;

use strict;
use warnings;

use Genome;
use File::Basename;

my $DEFAULT = '2.74';

class Genome::Model::Tools::Bsmap {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of BSMAP to use, default is $DEFAULT" },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run BSMAP or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools bsmap ...    
EOS
}

sub help_detail {                           
    return <<EOS 
More information about the BSMAP suite of tools can be found at http://code.google.com/p/bsmap/.
EOS
}

# TODO install this to $ENV{GENOME_SW}/bsmap/bsmap-2.1/
my %BSMAP_VERSIONS = (
    '2.1' => '/gsc/bin/bsmap',
    '2.43' => '/gscuser/cmiller/usr/src/bsmap-2.43/bsmap',
    '2.6' => '/gscuser/cmiller/usr/src/bsmap-2.6/bsmap',
    '2.6mod' => '/gscuser/cmiller/usr/src/bsmap-2.6.mod/bsmap',
    '2.7beta' => '/gscuser/cmiller/usr/src/bsmap-2.7beta/bsmap',
    '2.74' => '/gsc/bin/bsmap',
    'bsmap'   => 'bsmap',
);


sub bsmap_path {
    my $self = $_[0];
    return $self->path_for_bsmap_version($self->use_version);
}

sub available_bsmap_versions {
    my $self = shift;
    return keys %BSMAP_VERSIONS;
}

sub path_for_bsmap_version {
    my $class = shift;
    my $version = shift;
    unless (defined($version)) {
        $class->status_message("No version specified! Using default version '$DEFAULT'.");
        $version = $DEFAULT;
    }
    if (defined $BSMAP_VERSIONS{$version}) {
        return $BSMAP_VERSIONS{$version};
    }
    die('No path for bsmap version '. $version);
}

sub default_bsmap_version {
    die "default bsmap version: $DEFAULT is not valid" unless $BSMAP_VERSIONS{$DEFAULT};
    return $DEFAULT;
}
        
sub default_version { return default_bsmap_version; }


# TODO i'm guessing the following are not required

#sub supports_bam_input {
#    my $class = shift;
#    my $version = shift;
#
#    my %ok_versions = ();
#
#    return (exists $ok_versions{$version});
#
#}
#
#sub supports_multiple_reference {
#    my $class = shift;
#    my $version = shift;
#
#    my %ok_versions = ('0.5.9-pem0.1' => 1);
#
#    return exists $ok_versions{$version};
#}

1;

