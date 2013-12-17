package Genome::Model::Tools::Breakdancer;

use warnings;
use strict;

use Genome;

my $DEFAULT_VERSION = '2010_07_19';

class Genome::Model::Tools::Breakdancer{
    is => 'Command',
    has => [
        use_version => {
            is => 'Version',
            is_optional => 1,
            is_input => 1,
            default_value => $DEFAULT_VERSION,
            doc => "Version of breakdancer to use"
        },
     ],
};

my %BREAKDANCER_VERSIONS = (
    '0.0.1r59'   => {
        dir => '/gsc/scripts/pkg/bio/breakdancer/breakdancer-0.0.1r59',
        cfg => 'bam2cfg.pl',
        max => 'BreakDancerMax.pl',
    },
    '2010_02_17' => {
        dir => '/gsc/scripts/pkg/bio/breakdancer/breakdancer-2010_02_17/bin',
        cfg => 'bam2cfg.pl',
        max => 'BreakDancerMax.pl',
    },
    '2010_03_02' => {
        dir => '/gsc/scripts/pkg/bio/breakdancer/breakdancer-2010_03_02/bin',
        cfg => 'bam2cfg.pl',
        max => 'BreakDancerMax.pl',
    },
    '2010_06_24' => {
        dir => $ENV{GENOME_SW} . '/breakdancermax/breakdancer-20100624',
        cfg => 'perl/bam2cfg_2.pl',  #original bam2cfg.pl has problem close pipe filehandler
        max => 'cpp/breakdancer_max',
    },
    '2010_07_19' => {
        dir => $ENV{GENOME_SW} . '/breakdancermax/breakdancer-20100719',
        cfg => 'perl/bam2cfg_2.pl',
        max => 'cpp/breakdancer_max',
    },
    '1.0' => {
        dir => '/',
        cfg => 'usr/lib/breakdancer-max1.0/bam2cfg.pl',
        max => 'usr/bin/breakdancer-max1.0',
    },
    '1.1' => {
        dir => '/',
        cfg => 'usr/lib/breakdancer-max1.1/bam2cfg.pl',
        max => 'usr/bin/breakdancer-max1.1',
    },
    '1.2' => {
        dir => '/',
        cfg => 'usr/lib/breakdancer-max1.2/bam2cfg.pl',
        max => 'usr/bin/breakdancer-max1.2',
    },
    '1.3' => {
        dir => '/',
        cfg => 'usr/lib/breakdancer-max1.3/bam2cfg.pl',
        max => 'usr/bin/breakdancer-max1.3',
    },
    '1.3.7' => {
        dir => '/',
        cfg => 'usr/lib/breakdancer-max1.3.7/bam2cfg.pl',
        max => 'usr/bin/breakdancer-max1.3.7',
    },
    '1.4.1' => {
        dir => '/',
        cfg => 'usr/lib/breakdancer-max1.4.1/bam2cfg.pl',
        max => 'usr/bin/breakdancer-max1.4.1',
    },
    '1.4.2' => {
        dir => '/',
        cfg => 'usr/lib/breakdancer-max1.4.2/bam2cfg.pl',
        max => 'usr/bin/breakdancer-max1.4.2',
    },
);


sub help_brief {
    "discovers structural variation using breakdancer",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
EOS
}

sub help_detail {                           
    return <<EOS 
This tool discovers structural variation.  It generates an appropriate configuration based on
the input BAM files and then uses that configuration to run breakdancer.
EOS
}

sub breakdancer_path {
    my $self = shift;
    return $self->path_for_breakdancer_version($self->use_version);
}

sub breakdancer_max_command { 
    my $self = shift;
    return $self->breakdancer_max_command_for_version($self->use_version);
}

sub breakdancer_config_command { 
    my $self = shift;
    return $self->breakdancer_config_command_for_version($self->use_version);
}

sub available_breakdancer_versions {
    my $self = shift;
    return keys %BREAKDANCER_VERSIONS;
}

sub path_for_breakdancer_version {
    my ($self, $version) = @_;

    if (defined $BREAKDANCER_VERSIONS{$version}) {
        my $dir = $BREAKDANCER_VERSIONS{$version}->{dir};
        unless (-d $dir) {
            $self->error_message("breakdancer base dir $dir for version $version is not valid");
            die $self->error_message;
        }
        return $dir;
    }
    die 'No path for breakdancer version '. $version;
}

sub breakdancer_max_command_for_version {
    my ($self, $version) = @_;

    if (defined $BREAKDANCER_VERSIONS{$version}->{max}) {
        my $max_cmd = $self->path_for_breakdancer_version($version) . "/" .  $BREAKDANCER_VERSIONS{$version}->{max};
        unless (-s $max_cmd and -x $max_cmd) {
            $self->error_message("breakdancer_max command $max_cmd for version $version is not valid");
            die $self->error_message;
        }
        return $max_cmd;
    }
    die 'No breakdancer max command for breakdancer version '. $version;
}

sub breakdancer_config_command_for_version {
    my ($self, $version) = @_;

    if (defined $BREAKDANCER_VERSIONS{$version}->{cfg}) {
        my $cfg_cmd = $self->path_for_breakdancer_version($version) . "/" .  $BREAKDANCER_VERSIONS{$version}->{cfg};
        unless (-s $cfg_cmd and -x $cfg_cmd) {
            $self->error_message("breakdancer config command $cfg_cmd for version $version is not valid");
            die $self->error_message;
        }
        return join(' ', $^X, $cfg_cmd);
    }
    die 'No breakdancer config command for breakdancer version '. $version;
}

sub default_breakdancer_version {
    die "default breakdancer version: $DEFAULT_VERSION is not valid" unless $BREAKDANCER_VERSIONS{$DEFAULT_VERSION};
    return $DEFAULT_VERSION;
}
 
1;
