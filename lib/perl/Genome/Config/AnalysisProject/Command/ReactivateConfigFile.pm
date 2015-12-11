package Genome::Config::AnalysisProject::Command::ReactivateConfigFile;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::ReactivateConfigFile {
    is => 'Command::V2',
    has_input => [
        profile_item => {
            is => 'Genome::Config::Profile::Item',
            doc => 'Specified config file will be reenabled',
            shell_args_position => 1,
        }
    ],
    doc => "set the status of a profile item to 'active'"
};

sub execute {
    my $self = shift;
    $self->profile_item->status('active');
    return 1;
}

1;
