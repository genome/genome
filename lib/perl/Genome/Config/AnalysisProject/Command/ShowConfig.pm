package Genome::Config::AnalysisProject::Command::ShowConfig;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::ShowConfig {
    is => 'Genome::Object::Command::List',
    has_input => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'The AnalysisProject whose configuration to list',
            shell_args_position => 1,
        },
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::Config::Profile::Item',
        },
        show => {
            default_value => 'id,file_path,updated_at,is_concrete,analysis_menu_item.name,status,tags.name',
        },
    ],
    has_optional_input => [
        filter => {
            shell_args_position => 2,
        },
    ],
};

sub _base_filter {
    my $self = shift;

    return sprintf('analysis_project_id="%s"', $self->analysis_project->id);
}

1;
