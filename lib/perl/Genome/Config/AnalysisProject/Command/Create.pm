package Genome::Config::AnalysisProject::Command::Create;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::Create {
    is => 'Command::V2',
    has_input => [
        name => {
            is                  => 'Text',
            doc                 => 'the name of the analysis project to create',
            shell_args_position => 1
        },
        is_cle => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'If specified, this will flag this analysis project as being CLIA-related. All models will be created accordingly.',
        },
        is_production => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'If specified, this will flag this analysis project as being production data (and not CLIA-related). All models will be created accordingly.',
        },
        no_config => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'If specified, you will not be prompted for any menu items and the project will be created without config.'
        },
        analysis_menu_items => {
            is                  => 'Genome::Config::AnalysisMenu::Item',
            is_optional         => 1,
            is_many             => 1,
            doc                 => 'the analysis menu items to associate with this project.',
            require_user_verify => 1,
        },
    ],
};

sub help_brief {
    return 'Create a new analysis project';
}

sub help_synopsis {
    return "genome config analysis-project create 'my new project'";
}

sub help_detail {
    return <<"EOS"
Given a project name, this will set up a new analysis project.
It will create a new, empty, configuration set in the process.
EOS
}

sub execute {
    my $self = shift;

    if($self->is_cle and $self->is_production) {
        die('The --is-cle and --is-production options are mutually exclusive.');
    }

    $self->_validate_name($self->name);
    my @params = (
        name => $self->name,
        is_cle => $self->is_cle,
    );
    unless($self->is_production) {
        push @params, run_as => Genome::Sys->username;
    }
    my $project = Genome::Config::AnalysisProject->create(@params);
    die('No project created!') unless $project;

    eval {
        my $menu_items = $self->_get_menu_items();
        my $config_items = $self->_add_config_items_to_project($project, $menu_items);
    };
    if (my $error = $@) {
        $project->delete() if $project;
        $self->error_message('Failed to create Analysis Project!');
        die($error);
    }
    $self->status_message(q{Successfully created '} . $project->name . q{' (} . $project->id . q{).});
    return $project;
}

sub _validate_name {
    my $self = shift;
    my $name = shift;

    my $existing = Genome::Config::AnalysisProject->get(name => $name);
    if ($existing) {
        die "An Analysis Project with the name " . $name
        . " already exists! Please choose a unique name";
    }
}

sub _get_menu_items {
    my $self = shift;

    if ($self->analysis_menu_items()) {
        return [$self->analysis_menu_items];
    } elsif ($self->no_config()) {
        return [];
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
