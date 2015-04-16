package Genome::Config::AnalysisProject::Command::AddConfigFile;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::AddConfigFile {
    is => 'Genome::Config::AnalysisProject::Command::AddConfigBase',
    has_input => [
       config_file  => {
            is                  => 'Path',
            doc                 => 'path to the config file',
            shell_args_position => 2,
        },
        store_only => {
            is => 'Boolean',
            doc => 'If set, the config file will only be stored. (It will not be used for processing.)',
            default_value => 0,
        },
    ],
};

sub help_brief {
    return 'add a yml config file to an existing analysis project';
}

sub help_synopsis {
    return "genome config analysis-project add-config-file <analysis-project> <file-path>";
}

sub help_detail {
    return <<"EOS"
This command is used to link a custom configuration file directly to an analysis-project.

(See also `genome config analysis-project add-menu-item` for adding a standard configuration.)
EOS
}

sub valid_statuses {
    return ("Pending", "Hold", "In Progress", "Template");
}

sub _create_profile_items {
    my $self = shift;

    my $status = $self->store_only ? 'disabled' : 'active';

    my $result = Genome::Config::Profile::Item->create_from_file_path(
        file_path => $self->config_file,
        analysis_project => $self->analysis_project,
        status => $status,
    );
}

1;
