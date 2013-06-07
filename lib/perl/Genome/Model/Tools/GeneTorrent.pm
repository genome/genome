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
            default_value => '-W 72:00 -q lims-datatransfer',
        },
    ]
};

sub execute {
    my $self = shift;

    my $cmd = 'GeneTorrent'
        . ' --credential-file /gscuser/kochoa/mykey.pem'
        . ' --download https://cghub.ucsc.edu/cghub/data/analysis/download/' . $self->uuid
        . ' --path ' . $self->target_path
        . ' --log stdout:verbose'
        . ' --verbose 2'
        . ' --max-children 1'
        . ' --rate-limit 10';

    $self->status_message('Cmd: ' . $cmd);

    my $res = eval{ Genome::Sys->shellcmd(cmd => $cmd); };

    if ( not $res ) {
        $self->error_message('Cannot execute command "' . $cmd . '" : ' . $@);
        return;
    }

    return 1;
}

1;
