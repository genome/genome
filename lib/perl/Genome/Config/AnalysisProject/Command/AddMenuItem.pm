package Genome::Config::AnalysisProject::Command::AddMenuItem;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::AddMenuItem {
    is => 'Command::V2',
    has_input => [
       analysis_project  => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis project to add the config file to',
            shell_args_position => 1,
        },
        analysis_menu_items => {
            is                  => 'Genome::Config::AnalysisMenu::Item',
            is_optional         => 1,
            is_many             => 1,
            doc                 => 'the analysis menu items to associate with this project.',
            require_user_verify => 1,
        }
    ],
};

sub help_brief {
    return 'add a menu item to an existing analysis project';
}

sub help_synopsis {
    return "genome config analysis-project add-config-config-file <analysis-project> <file-path>";
}

sub help_detail {
    return <<"EOS"
Given an analysis project and a config file, this will associate the two
EOS
}

sub execute {
    my $self = shift;

    eval {
        my $menu_items = $self->_get_menu_items();
        my $config_items = $self->_add_config_items_to_project($self->analysis_project, $menu_items);
    };
    if (my $error = $@) {
        $self->error_message('Failed to add to Analysis Project!');
        die($error);
    }

    return $self->analysis_project;
}

sub _get_menu_items {
    my $self = shift;

    if ($self->analysis_menu_items()) {
        return [$self->analysis_menu_items];
    } else {
        my $class_name = 'Genome::Config::AnalysisMenu::Item';
        return [$self->resolve_param_value_from_cmdline_text(
                {
                    name => 'analysis_menu_items',
                    class => $class_name,
                    value => [$class_name->get()],
                }
            )];
    }
}

sub _add_config_items_to_project {
    my $self = shift;
    my $project = shift;
    my $menu_items = shift;

    for (@$menu_items) {
        Genome::Config::Profile::Item->create(
            analysis_menu_item => $_,
            analysis_project => $project,
        );
    }

    return 1;
}

1;
