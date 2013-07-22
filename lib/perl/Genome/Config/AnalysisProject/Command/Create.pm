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
        analysis_menu_item => {
            is          => 'Genome::Config::AnalysisMenuItem',
            is_optional => 1,
            doc         => 'the analysis menu item to associate with this project. default value is all available processing',
            default     => '702fb86975c64ffb8b965898e7c4c5ad'
        }
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

    my $existing = Genome::Config::AnalysisProject->get(name => $self->name);
    if ($existing) {
        die "An Analysis Project with the name " . $self->name
            . " already exists! Please choose a unique name";
    }

    my $ap;
    eval {
        $ap = Genome::Config::AnalysisProject->create(
            name => $self->name,
            _analysis_menu_item_id => $self->analysis_menu_item->id,
        );
    };
    if (my $error = $@) {
        $ap->delete();
        $self->error_message('Failed to create Analysis Project!');
        die($error);
    }

    UR::Context->commit();
    $self->status_message('Successfully created analysis project ' . $ap->name .
        ' ID: ' . $ap->id);

    return $ap;
}

1;
