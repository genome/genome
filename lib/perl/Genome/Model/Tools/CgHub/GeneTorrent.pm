package Genome::Model::Tools::CgHub::GeneTorrent;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::CgHub::GeneTorrent {
    is => 'Genome::Model::Tools::CgHub::Base',
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
    ],
    has_calculated => {
        source_url => {
            calculate_from => [qw/ uuid /],
            calculate => q( return 'https://cghub.ucsc.edu/cghub/data/analysis/download/'.$uuid; ),
        },
        credential_file => { #FIXME each user should have their own
            calculate => q( return '/gscuser/kochoa/mykey.pem'; ),
        }
    },
    doc => 'Download files from CG Hub using gene-torrent',
};

sub _build_command {
    my $self = shift;

    my $cmd = 'gtdownload'
        . ' --credential-file '.$self->credential_file 
        . ' --download ' . $self->source_url
        . ' --path ' . $self->target_path
        . ' --log stdout:verbose'
        . ' --verbose 2'
        . ' --max-children 2'
        . ' --rate-limit '.$self->rate_limit # mega-BYTES per second (see internet_download_mbps above)
        . ' --inactivity-timeout ' . 3 * 60 * 24   # in minutes - instead of bsub -W
    ;
}

sub _verify_success {
    my $self = shift;

    if ( not -s $self->target_path ) {
        $self->error_message('Failed to download source URL! '.$self->source_url);
        return;
    }

    return 1;
}

sub rate_limit {
    my $self = shift;
    $self->lsf_resource =~ /internet_download_mbps=(\d+)/;
    my $internet_download_mbps = $1;
    return 10 if not $internet_download_mbps; # previously hard coded
    return int($internet_download_mbps / 8); # convert to bytes
}

1;

