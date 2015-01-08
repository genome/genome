package Genome::Config::AnalysisProject::Command::AddConfigFile;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::AddConfigFile {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
       config_file  => {
            is                  => 'Path',
            doc                 => 'path to the config file',
            shell_args_position => 2,
        },
        store_only => {
            is => 'Boolean',
            doc => 'If set, the config file will only be stored (it will not be used for processing).  Defaults to 0',
            default_value => 0,
        },
        reprocess_existing => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Reprocess any existing instrument data with the new config',
        },
    ],
    has_optional_input => [
        tag => {
            is => 'Genome::Config::Tag',
            is_many => 1,
            doc => 'Optional tags to add',
        }
    ]
};

sub help_brief {
    return 'add a yml config file to an existing analysis project';
}

sub help_synopsis {
    return "genome config analysis-project add-config-config-file <analysis-project> <file-path>";
}

sub help_detail {
    return <<"EOS"
Given an analysis project and a config file, this will associate the two
EOS
}

sub valid_statuses {
    return ("Pending", "Hold", "In Progress", "Template");
}

sub execute {
    my $self = shift;
    my $status = $self->store_only ? 'disabled' : 'active';

    my $result = Genome::Config::Profile::Item->create_from_file_path(
        file_path => $self->config_file,
        analysis_project => $self->analysis_project,
        status => $status,
    );

    $self->_apply_tags($result);

    if($self->reprocess_existing){
        $self->_mark_instrument_data_bridges;
    }

    return $result;
}

sub _mark_instrument_data_bridges {
    my $self = shift;
    my $analysis_project = $self->analysis_project;
    for my $bridge ($analysis_project->analysis_project_bridges){
        $bridge->status('new');
    }
    return 1;
}

sub _apply_tags {
    my ($self, $profile_item) = @_;
    for my $tag ($self->tag){
        $profile_item->add_tag($tag);
    }
    return 1;
}

1;
