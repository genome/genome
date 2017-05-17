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

    my %inputs;

    for my $input ($build->inputs) {
        my $key = $input->name;
        $inputs{$key} = $build->$key;
        #TODO handle complex inputs
    }

    my $yaml = File::Spec->join($build->data_directory, 'inputs.yaml');
    YAML::DumpFile($yaml, \%inputs);
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
    my $jobstore_dir = File::Spec->join($build->data_directory, 'jobstore');
    my $log_file = File::Spec->join($build->log_directory, 'toil.log');

    Genome::Sys->shellcmd(
        cmd => [
            'cwltoil',
            '--disableCaching', '--stats', '--logLevel=DEBUG',
            '--workDir', $tmp_dir,
            '--jobStore', $jobstore_dir,
            "--logFile=$log_file",
            "--outdir=$results_dir",
            '--batchSystem', 'lsf',
            $build->main_workflow_file,
            $yaml
        ],
        input_files => [$build->main_workflow_file, $yaml],
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
