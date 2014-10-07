package Genome::Model::Tools::CgHub::GeneTorrent;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::CgHub::GeneTorrent {
    is => "Genome::Model::Tools::CgHub::Base",
    has_input => [
        uuid => {
            is => "Text",
            is_output => 1,
        },
        target_path => {
            is => "Text",
            is_output => 1,
        },
    ],
    has => [
        lsf_resource => {
            # mbps -> mega-BITS per second (see --rate-limit below)
            default_value => '-q lims-long -R "rusage[internet_download_mbps=80]"',
        },
    ]
};

sub _build_command {
    my $self = shift;

    # version 3.3.4 has GeneTorrent binary
    # version 3.8.3 has gtdownload binary
    my $run_gtdownload = eval{ 
        local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 0;
        $self->_run_command("gtdownload --help > /dev/null");
    };
    my $exe = ( $run_gtdownload ) ? 'gtdownload' : 'GeneTorrent';

    return "$exe"
        . ' --credential-file /gscuser/kochoa/mykey.pem'    # TODO: do not hardcode
        . ' --download https://cghub.ucsc.edu/cghub/data/analysis/download/' . $self->uuid
        . ' --path ' . $self->target_path
        . ' --log stdout:verbose'
        . ' --verbose 2'
        . ' --max-children 2'
        . ' --rate-limit '.$self->rate_limit # mega-BYTES per second (see internet_download_mbps above)
        . ' --inactivity-timeout ' . 3 * 60 * 24   # in minutes - instead of bsub -W
    ;
}

sub _verify_success { return 1; }

sub rate_limit {
    my $self = shift;
    $self->lsf_resource =~ /internet_download_mbps=(\d+)/;
    my $internet_download_mbps = $1;
    return 10 if not $internet_download_mbps; # previously hard coded
    return int($internet_download_mbps / 8); # convert to bytes
}

1;

