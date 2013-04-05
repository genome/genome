
package Genome::Model::Tools::Tophat;

use strict;
use warnings;

use Genome;
use File::Basename;
use IPC::Cmd;

my $DEFAULT = '1.4.0';

class Genome::Model::Tools::Tophat {
    is => 'Command',
    has => [
        use_version => { is => 'Version', is_optional => 1, default_value => $DEFAULT, doc => "Version of tophat to use, default is $DEFAULT" },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to run Tophat or work with its output files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools tophat ...    
EOS
}

sub help_detail {
    return <<EOS
More information about the Tophat aligner can be found at http://tophat.cbcb.umd.edu/.
EOS
}


my %TOPHAT_VERSIONS = (
    '0.7.1'  => $ENV{GENOME_SW} . '/tophat/tophat-0.7.1-64/bin/tophat',
    '0.7.2'  => $ENV{GENOME_SW} . '/tophat/tophat-0.7.2-64/bin/tophat',
    '1.0.10' => $ENV{GENOME_SW} . '/tophat/tophat-1.0.10-64/bin/tophat',
    # These are 64-bit installations.  tophat is really a python script and is OS dependent.
    '1.0.12' => $ENV{GENOME_SW} . '/tophat/tophat-1.0.12/bin/tophat',
    '1.0.13' => $ENV{GENOME_SW} . '/tophat/tophat-1.0.13/bin/tophat',
    '1.0.14' => $ENV{GENOME_SW} . '/tophat/tophat-1.0.14/tophat',
    '1.1.0'  => $ENV{GENOME_SW} . '/tophat/tophat-1.1.0/tophat',
    '1.1.2'  => $ENV{GENOME_SW} . '/tophat/tophat-1.1.2/tophat',
    '1.1.4'  => $ENV{GENOME_SW} . '/tophat/tophat-1.1.4/tophat',
    '1.2.0'  => $ENV{GENOME_SW} . '/tophat/tophat-1.2.0/tophat',
    '1.3.0'  => $ENV{GENOME_SW} . '/tophat/tophat-1.3.0/tophat',
    '1.3.1'  => $ENV{GENOME_SW} . '/tophat/tophat-1.3.1/tophat',
    'tophat' => 'tophat',
);


sub tophat_path {
    my $self = $_[0];
    return $self->path_for_tophat_version($self->use_version);
}

sub available_tophat_versions {
    my $self = shift;
    my @legacy_versions = keys %TOPHAT_VERSIONS;
    my @local_versions = qx/ update_alternatives --list tophat /;
    @local_versions = map { /^.*tophat(.*)/; $1; } @local_versions;
    return @legacy_versions, @local_versions;
}

sub path_for_tophat_version {
    my $class = shift;
    my $version = shift;

    my $path = IPC::Cmd::can_run("tophat" . $version);
    return $path if ($path);

    if (defined $TOPHAT_VERSIONS{$version}) {
        return $TOPHAT_VERSIONS{$version};
    }
    die('No path for tophat version '. $version);
}

sub default_tophat_version {
    die "default tophat version: $DEFAULT is not valid" unless __PACKAGE__->path_for_tophat_version($DEFAULT);
    return $DEFAULT;
}

1;

