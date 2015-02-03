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
    has_optional_input => [
        tags => {
            is                  => 'Boolean',
            doc                 => 'copy tags (if any) that are attached to the configs',
            default_value       => 0,
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

    my @configs_to_copy = grep { $_->status ne 'disabled' } $from_project->config_items;
    $self->_copy_config_profile_items_to_project($to_project, @configs_to_copy);

    return 1;
}

sub _copy_config_profile_items_to_project {
    my $self = shift;
    my $analysis_project = shift;
    my @config_profile_items = @_;

    for my $original_config_item (@config_profile_items) {
        my $new_item;

        my @params = (
            analysis_project => $analysis_project,
            status => $original_config_item->status,
        );
        if($self->tags) {
            push @params, tags => [$original_config_item->tags];
        }

        if ($original_config_item->analysis_menu_item) {
            $new_item = Genome::Config::Profile::Item->create(
                analysis_menu_item => $original_config_item->analysis_menu_item,
                @params,
            );
        } else {
            $new_item = Genome::Config::Profile::Item->create_from_file_path(
                file_path => $original_config_item->file_path,
                @params,
            );
        }
    }

    return 1;
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__(@_);
    my $status = $self->to_project->status;
    unless(grep{$_ eq $status} ("Pending", "Hold", "In Progress", "Template")){
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['to_project'],
            desc => "Can't add config file to analysis project with status: $status"
        );
    }
    return @errors;
}

1;
