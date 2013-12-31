package Genome::Config::AnalysisProject::Command::CopyConfig;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::CopyConfig {
    is => 'Command::V2',
    has_input => [
        from_project => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'source project of config to copy',
            shell_args_position => 1
        },
        to_project => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'destination project of config to copy',
            shell_args_position => 2
        },
    ],
};

sub help_brief {
    return 'Copies the configuration profile from one analysis project to another';
}

sub help_synopsis {
    return "genome config analysis-project copy-config FROM_PROJECT TO_PROJECT";
}

sub help_detail {
    return <<"EOS"
This command will copy the configuration profile items from one analysis project to another.
If the project has menu items selected, it will use those menu items in "non-concrete" form.
Any custom configuration files will immediately be copied over in concrete form.
EOS
}

sub execute {
    my $self = shift;

    my $from_project = $self->from_project;
    my $to_project = $self->to_project;

    $self->_copy_config_profile_items_to_project($to_project, $from_project->config_items);

    return 1;
}

sub _copy_config_profile_items_to_project {
    my $self = shift;
    my $analysis_project = shift;
    my @config_profile_items = @_;

    for my $original_config_item (@config_profile_items) {
        if ($original_config_item->analysis_menu_item) {
            Genome::Config::Profile::Item->create(
                analysis_project => $analysis_project,
                analysis_menu_item => $original_config_item->analysis_menu_item,
            );
        } else {
            Genome::Config::Profile::Item->create_from_file_path(
                analysis_project => $analysis_project,
                file_path => $original_config_item->file_path,
            );
        }
    }

    return 1;
}

1;
