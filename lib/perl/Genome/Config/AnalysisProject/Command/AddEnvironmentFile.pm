package Genome::Config::AnalysisProject::Command::AddEnvironmentFile;

use strict;
use warnings;

use File::Spec;

use Genome;

class Genome::Config::AnalysisProject::Command::AddEnvironmentFile {
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
    return 'add a yml environment file to an existing analysis project';
}

sub help_synopsis {
    return "genome config analysis-project add-environment-file <analysis-project> <file-path>\n";
}

sub help_detail {
    return <<"EOS"
This command is used to link an environment configuration file to an analysis-project.  This specifies properties of the environment in which builds for the analysis-project will run.
EOS
}

sub valid_statuses {
    return ('Pending', 'In Progress', 'Hold');
}

sub execute {
    my $self = shift;
    my $analysis_project = $self->analysis_project;
    if($analysis_project->environment_config_dir) {
        $self->fatal_message('Project %s already has an existing environment file.', $analysis_project->__display_name__);
    }

    my $allocation = $analysis_project->disk_allocations;
    unless ($allocation) {
        $allocation = Genome::Disk::Allocation->create(
            owner_id => $analysis_project->id,
            owner_class_name => $analysis_project->class,
            disk_group_name => Genome::Config::get('disk_group_references'),
            allocation_path => 'analysis_project/' . $analysis_project->id,
            kilobytes_requested => ((-s $self->environment_file)/1024 + 1),
        );
    }

    if($ENV{UR_DBI_NO_COMMIT}) {
        my $config_subpath = Genome::Config::config_subpath;
        my $filename = $config_subpath->basename;

        my $dir = File::Spec->join($allocation->absolute_path, $config_subpath->dir->stringify);
        Genome::Sys->create_directory($dir);

        my $installed_path = File::Spec->join($dir, $filename);
        Genome::Sys->copy_file($self->environment_file, $installed_path);
        $allocation->reallocate();

        $self->status_message('Environment file installed to %s', $installed_path);
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

        $self->status_message('Environment file queued for installation');
    }
    return 1;
}

1;
