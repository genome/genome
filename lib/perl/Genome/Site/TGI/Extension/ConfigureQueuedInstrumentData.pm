package Genome::Config::Command::ConfigureQueuedInstrumentData;

use strict;
use warnings;

use Genome;
use Genome::Config::Command::ConfigureQueuedInstrumentData; #load real module first

use Sub::Install;

Sub::Install::reinstall_sub ({
    as => '_import_instrument_data',
    code => \&_import_lims_data,

});

sub _import_lims_data {
    my $self = shift;
    my $current_pair = shift;

    my @cmd = ($^X);
    push @cmd, '-I', $_ for UR::Util::used_libs;
    push @cmd,
        '-MGenome',
        '-e', 'Genome::Site::TGI::Command->execute_with_shell_params_and_exit()',
        'import-data-from-lims',
        '--instrument-data', $current_pair->instrument_data_id,
        '--analysis-project', $current_pair->analysis_project_id,
    ;

    Genome::Sys->shellcmd(cmd => \@cmd);
}

1;
