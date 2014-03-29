package Genome::Config::AnalysisProject::Command::DisableConfigFile;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::DisableConfigFile {
    is => 'Command::V2',
    has_input => [
        profile_item => {
            is => 'Genome::Config::Profile::Item',
            doc => 'Specified config file will no longer be used',
            shell_args_position => 1,
        }
    ],
};

sub execute {
    my $self = shift;
    $self->profile_item->status('disabled');
    return 1;
}

1;
