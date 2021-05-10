package Genome::Config::AnalysisProject::Command::UpdateEnvironmentFile;

use strict;
use warnings;

use File::Spec;

use Genome;

class Genome::Config::AnalysisProject::Command::UpdateEnvironmentFile {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
        environment_file  => {
            is                  => 'Path',
            doc                 => 'path to the environment config file',
            shell_args_position => 2,
        },
    ],
};

sub help_brief {
    return "update an analysis project's existing yml environment file";
}

sub help_synopsis {
    return "genome config analysis-project update-environment-file <analysis-project> <file-path>\n";
}

sub help_detail {
    return <<"EOS"
This command is used to update the existing environment configuration file for an analysis-project.
EOS
}

sub valid_statuses {
    return ('Pending', 'In Progress', 'Hold');
}

sub execute {
    my $self = shift;
    my $analysis_project = $self->analysis_project;
    unless($analysis_project->environment_config_dir) {
        $self->fatal_message('Project %s has no existing environment file.', $analysis_project->__display_name__);
    }

    my $allocation = $analysis_project->disk_allocations;

    if($ENV{UR_DBI_NO_COMMIT}) {
        my $config_subpath = Genome::Config::config_subpath;
        my $filename = $config_subpath->basename;

        my $dir = File::Spec->join($allocation->absolute_path, $config_subpath->dir->stringify);

        my $installed_path = File::Spec->join($dir, $filename);
        Genome::Sys->copy_file($self->environment_file, $installed_path);
        $allocation->reallocate();

        $self->status_message('Environment file updated at %s', $installed_path);
    }
    else {
        my $env_filename = 'env.' . $analysis_project->id . '.yaml';
        Genome::Sys->shellcmd(
            cmd => [
                '/usr/bin/gsutil/gsutil',
                'cp',
                $self->environment_file,
                'gs://gms_environment_config/' . $env_filename,
            ],
            input_files => [$self->environment_file],
        );

        $self->status_message('Environment file queued for update');
    }
    return 1;
}

1;
