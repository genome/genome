package Genome::Model::CwlPipeline::Command::Run;

use strict;
use warnings;

use File::Spec;
use Genome;
use Genome::Utility::File::Mode qw();
use YAML;

class Genome::Model::CwlPipeline::Command::Run {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build::CwlPipeline',
        },
    ],
    doc => 'wrapper command to run "cwltoil"'
};

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self = shift;

    my $yaml = $self->prepare_yaml;
    my ($tmp_dir, $results_dir) = $self->prepare_directories;
    $self->run_toil($yaml, $tmp_dir, $results_dir);
    $self->preserve_results($results_dir);
    $self->cleanup($tmp_dir);

    return 1;
}

sub prepare_yaml {
    my $self = shift;
    my $build = $self->build;

    my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');

    my $pp = $build->processing_profile;
    my $da = Genome::Disk::Allocation->get(owner => $pp);
    my $dir = $da->absolute_path;
    my $preprocessing_script = File::Spec->join($dir, 'process_inputs.pl');

    if(-e $preprocessing_script) {
        Genome::Sys->shellcmd(
            cmd => [$^X, $preprocessing_script, $build->id],
            output_files => [$yaml],
        );
    } else {
        my %inputs;

        for my $input ($build->inputs) {
            my $key = $input->name;
            if ($input->value_class_name->isa('UR::Value')) {
                $inputs{$key} = $input->value_id;
            } else {
                #TODO handle complex inputs?
                $self->fatal_message('Cannot handle object inputs yet!');
            }
        }

        YAML::DumpFile($yaml, \%inputs);
    }

    return $yaml;
}

sub prepare_directories {
    my $self = shift;
    my $build = $self->build;

    my $data_directory = $build->data_directory;

    my $tmp_dir = File::Spec->join($data_directory, 'tmp');
    my $results_dir = File::Spec->join($data_directory, 'results');

    Genome::Sys->create_directory($tmp_dir);
    Genome::Sys->create_directory($results_dir);

    return ($tmp_dir, $results_dir);
}

sub run_toil {
    my $self = shift;
    my $yaml = shift;
    my $tmp_dir = shift;
    my $results_dir = shift;

    my $build = $self->build;
    my $model = $build->model;
    my $jobstore_dir = File::Spec->join($build->data_directory, 'jobstore');
    my $log_file = File::Spec->join($build->log_directory, 'toil.log');

    my @restart;
    if (-e $jobstore_dir) {
        push @restart, '--restart';
    }

    my $primary_docker_image = $model->primary_docker_image;
    local $ENV{LSB_SUB_ADDITIONAL} = $primary_docker_image;

    my $default_queue = Genome::Config::get('lsf_queue_build_worker_alt');
    local $ENV{LSB_DEFAULTQUEUE} = $default_queue;

    #toil relies on reading the output from bsub
    delete local $ENV{BSUB_QUIET};

    Genome::Sys->shellcmd(
        cmd => [
            'cwltoil',
            '--disableCaching', '--logLevel=DEBUG',
            @restart,
            '--workDir', $tmp_dir,
            '--jobStore', $jobstore_dir,
            "--logFile=$log_file",
            "--outdir=$results_dir",
            '--batchSystem', 'lsf',
            $model->main_workflow_file,
            $yaml
        ],
        input_files => [$model->main_workflow_file, $yaml],
    );
}

sub preserve_results {
    my $self = shift;
    my $results_dir = shift;

    for my $file ($results_dir, glob("$results_dir/*")) {
        Genome::Utility::File::Mode::mode($file)->rm_all_writable;
    }

    return 1;
}

sub cleanup {
    my $self = shift;
    my $tmp_dir = shift;

    return Genome::Sys->remove_directory_tree($tmp_dir);
}

1;
