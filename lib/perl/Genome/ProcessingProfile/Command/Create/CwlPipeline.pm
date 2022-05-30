package Genome::ProcessingProfile::Command::Create::CwlPipeline;

use strict;
use warnings;

use File::Spec;
use Genome;

class Genome::ProcessingProfile::Command::Create::CwlPipeline {
    is => 'Genome::ProcessingProfile::Command::Create::Base',
    has_input => [
        cwl_directory => {
            is => 'Directory',
            doc => 'a directory of CWL for the pipeline',
        },
        main_workflow_file => {
            is => 'Text',
            doc => 'name of the main workflow file within the directory',
        },
        primary_docker_image => {
            is => 'Text',
            doc => 'docker image for the main toil worker jobs',
            example_values => ['docker(mgibio/rnaseq)'],
        },
        short_pipeline_name => {
            is => 'Text',
            is_optional => 1,
            doc => 'short name for pipeline to include in default model names',
        },
    ],
    has_transient_optional => [
        based_on => {
            value => undef,
            doc => 'the based_on option is unavailable for this type of processing profile',
        },
    ],
};

sub _target_class_name {
    return 'Genome::ProcessingProfile::CwlPipeline';
}

sub execute {
    my $self = shift;

    $self->_verify_entrypoint_exists;

    my $ok = $self->SUPER::_execute_body;
    my $new_pp = $self->created_processing_profile;
    unless ($ok and $new_pp) {
        $self->fatal_message('Not uploading workflow files due to error creating processing profile.');
    }

    my $cwl_directory = $self->cwl_directory;
    my $da = Genome::Disk::Allocation->create(
        owner_id => $new_pp->id,
        owner_class_name => $new_pp->class,
        kilobytes_requested => Genome::Sys->disk_usage_for_path($cwl_directory),
        disk_group_name => Genome::Config::get('disk_group_references'),
        allocation_path => 'processing-profile/cwl-pipeline/' . $new_pp->id,
        skip_allocation_path_creation => 1,
    );

    if($ENV{UR_DBI_NO_COMMIT}) {
        #FIXME deallocate on rsync errors or on PP delete
        Genome::Sys->rsync_directory(
            source_directory => $cwl_directory,
            target_directory => $da->absolute_path,
        );
    }
    else {
        my $pp_dirname = 'pp.' . $new_pp->id;
        Genome::Sys->shellcmd(
            cmd => [
                '/usr/bin/python3',
                '/usr/bin/gsutil/gsutil',
                '-m', '-q', 'cp', '-r',
                $cwl_directory,
                Genome::Config::get('gcp_config_bucket') . $pp_dirname,
            ],
            input_directories => [$cwl_directory],
        );

        $self->status_message('Processing profile queued for installation');
    }

    $new_pp->main_workflow_file(
        File::Spec->join($da->absolute_path, $self->main_workflow_file)
    );

    return 1;
}

sub _verify_entrypoint_exists {
    my $self = shift;

    my $entrypoint = File::Spec->join($self->cwl_directory, $self->main_workflow_file);
    unless (-e $entrypoint) {
        $self->fatal_message('Specified main workflow file not found in CWL directory.  Does <%s> exist?', $entrypoint);
    }

    return 1;
}

1;
