package Genome::Model::Tools::GeneTorrent;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::GeneTorrent {
    is => "Command::V2",
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
#            default_value => '-W 72:00 -q lims-datatransfer -R "rusage[internet_download_mbps=80]"',
            # -W bsub option silently ignored
            # mbps -> megabits per second (see --rate-limit below)
            default_value => '-q seqmgr-long -R "rusage[internet_download_mbps=80]"',
        },
    ]
};

sub execute {
    my $self = shift;

    # genetorrent-download debian package for 10.04 lives in /cghub
    local $ENV{'PATH'} = $ENV{'PATH'} . ':/cghub/bin';

    # version 3.3.4 has GeneTorrent binary
    # version 3.8.3 has gtdownload binary
    my $exe = do {
        `gtdownload --help`;
        ($? == 0) ? 'gtdownload' : 'GeneTorrent';
    };
    my $cmd = "$exe"
        . ' --credential-file /gscuser/kochoa/mykey.pem'    # TODO: do not hardcode
        . ' --download https://cghub.ucsc.edu/cghub/data/analysis/download/' . $self->uuid
        . ' --path ' . $self->target_path
        . ' --log stdout:verbose'
        . ' --verbose 2'
        . ' --max-children 2'
        . ' --rate-limit 10'    # megabytes per second (see internet_download_mbps above)
        . ' --inactivity-timeout ' . 3 * 60 * 24   # in minutes - instead of bsub -W
    ;

    $self->debug_message('Cmd: ' . $cmd);

    my $res = eval{ Genome::Sys->shellcmd(cmd => $cmd); };

    if ( not $res ) {
        $self->error_message('Cannot execute command "' . $cmd . '" : ' . $@);
        return;
    }

    return 1;
}

1;
