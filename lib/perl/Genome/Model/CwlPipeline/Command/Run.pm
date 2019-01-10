package Genome::Model::CwlPipeline::Command::Run;

use strict;
use warnings;

use File::Spec;
use Genome;
use Genome::Utility::File::Mode qw();
use YAML;
use JSON qw(to_json);
use File::Compare qw();

class Genome::Model::CwlPipeline::Command::Run {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build::CwlPipeline',
        },
        lsf_resource => {
            is_param => 1,
            value => Genome::Config::get('lsf_resource_cwl_runner'),
        },
    ],
    doc => 'wrapper command to run "cwltoil"'
};

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self = shift;

    my $yaml = $self->prepare_yaml;
    my ($tmp_dir, $results_dir) = $self->prepare_directories;
    my $cwl_runner = Genome::Config::get('cwl_runner');
    if ($cwl_runner eq 'cromwell') {
        $self->run_cromwell($yaml, $tmp_dir, $results_dir);
        $self->cleanup($tmp_dir);
    } elsif ($cwl_runner eq 'toil') {
        $self->run_toil($yaml, $tmp_dir, $results_dir);
        $self->cleanup($tmp_dir, $results_dir);
    } else {
        $self->fatal_message('Unknown CWL runner: %s', $cwl_runner);
    }
    $self->preserve_results($results_dir);

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

sub run_cromwell {
    my $self = shift;
    my $yaml = shift;
    my $tmp_dir = shift;
    my $results_dir = shift;

    my $build = $self->build;
    my $model = $build->model;

    #cromwell relies on reading the output from bsub
    delete local $ENV{BSUB_QUIET};

    my $config_file = $self->_generate_cromwell_config($tmp_dir, $results_dir);
    my $labels_file = $self->_generate_cromwell_labels;

    my $truststore_file = Genome::Config::get('cromwell_truststore_file');
    my $truststore_auth = Genome::Config::get('cromwell_truststore_auth');

    Genome::Sys->shellcmd(
        cmd => [
            '/usr/bin/java',
            sprintf('-Dconfig.file=%s', $config_file),
            sprintf('-Djavax.net.ssl.trustStorePassword=%s', $truststore_auth),
            sprintf('-Djavax.net.ssl.trustStore=%s', $truststore_file),
            '-jar', '/opt/cromwell.jar',
            'run',
            '-t', 'cwl',
            '-l', $labels_file,
            '-i', $yaml,
            $model->main_workflow_file,
        ],
        input_files => [$model->main_workflow_file, $yaml],
    );

    $self->_stage_cromwell_outputs($results_dir);
}

sub _generate_cromwell_config {
    my $self = shift;
    my $tmp_dir = shift;

    my $build = $self->build;
    my $log_dir = $build->log_directory;

    my $primary_docker_image = $build->model->primary_docker_image;

    my $default_queue = Genome::Config::get('lsf_queue_build_worker_alt');

    my $server = Genome::Config::get('cromwell_server');
    my $auth = Genome::Config::get('cromwell_auth');
    my $user = Genome::Config::get('cromwell_user');

    my $config_file = File::Spec->join($build->data_directory, 'cromwell.config');

    my $config = <<'EOCONFIG'
include required(classpath("application"))

backend {
  default = "LSF"
  providers {
    LSF {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
        Int cpu = 1
        Int? memory_kb
        Int? memory_mb
        String? docker
        """

        submit = """
        bsub \
        -J ${job_name} \
        -cwd ${cwd} \
EOCONFIG
    ;
    $config .= <<EOCONFIG
        -o /dev/null \\
        -e $log_dir/cromwell-%J.err \\
        -a '$primary_docker_image' \\
        -q '$default_queue' \\
EOCONFIG
        ;
    $config .= <<'EOCONFIG'
        ${"-M " + memory_kb} \
        ${"-n " + cpu} \
        ${"-R \"select[mem>" + memory_mb + "] rusage[mem=" + memory_mb + "]\""} \
        /bin/bash ${script}
        """

        submit-docker = """
        LSF_DOCKER_VOLUMES=${cwd}:${docker_cwd} \
        bsub \
        -J ${job_name} \
        -cwd ${cwd} \
EOCONFIG
    ;
    $config .= <<EOCONFIG
        -o /dev/null \\
        -e $log_dir/cromwell-%J.err \\
        -q '$default_queue' \\
EOCONFIG
        ;
    $config .= <<'EOCONFIG'
        ${"-a \"docker(" + docker + ")\""} \
        ${"-M " + memory_kb} \
        ${"-n " + cpu} \
        ${"-R \"select[mem>" + memory_mb + "] rusage[mem=" + memory_mb + "]\""} \
        /bin/bash ${script}
        """

        kill = "bkill ${job_id}"
        docker-kill = "bkill ${job_id}"
        check-alive = "bjobs -noheader -o stat ${job_id} | /bin/grep 'PEND\\|RUN'"
        job-id-regex = "Job <(\\d+)>.*"
EOCONFIG
;
    $config .= <<EOCONFIG
        root = "$tmp_dir/cromwell-executions"
      }
    }
  }
}
workflow-options {
  workflow-log-dir = "$log_dir/cromwell-workflow-logs"
}
database {
  profile = "slick.jdbc.MySQLProfile\$"
  db {
    driver = "com.mysql.jdbc.Driver"
    url = "jdbc:$server"
    user = "$user"
    password = "$auth"
    connectionTimeout = 5000
    numThreads = 5
  }
}
EOCONFIG
;

    Genome::Sys->write_file($config_file, $config);
    return $config_file;
}

sub _generate_cromwell_labels {
    my $self = shift;
    my $build = $self->build;

    my $labels_file = File::Spec->join($build->data_directory, 'cromwell.labels');

    my $data = {
        build => $build->id,
        model => $build->model->id,
        analysis_project => $build->model->analysis_project->id,
    };

    Genome::Sys->write_file($labels_file, to_json($data));
    return $labels_file;
}

sub _stage_cromwell_outputs {
    my $self = shift;
    my $results_dir = shift;

    my $build = $self->build;

    my $results = Genome::Cromwell->query( [{ label => 'build:' . $build->id }] );
    if ($results->{totalResultsCount} != 1) {
        $self->fatal_message('Failed to find workflow.  Got: %s', $results);
    }

    my $workflow_id = $results->{results}->[0]->{id};
    my $output_result = Genome::Cromwell->outputs($workflow_id);

    my $outputs = $output_result->{outputs};

    for my $output_name (keys %$outputs) {
        my $info = $outputs->{$output_name};
        
        $self->_stage_cromwell_output($results_dir, $info); 
    }

    return 1;
}

sub _stage_cromwell_output {
    my $self = shift;
    my $results_dir = shift;
    my $item = shift;

    if (ref $item eq 'ARRAY') {
        map $self->_stage_cromwell_output($results_dir, $_), @$item;
        return 1;
    }

    return unless ref $item eq 'HASH';
    my $location = $item->{location};
    return unless $location;

    my @secondary = map { $_->{location} } @{ $item->{secondaryFiles} };

    for my $source ($location, @secondary) {
        my (undef, $dir, $file) = File::Spec->splitpath($source);

        my $destination = File::Spec->join($results_dir, $file);
        if (-e $destination) {
            if( File::Compare::compare($source, $destination) == 0 ) {
                $self->warning_message('Skipping staging of %s--an identical file was already staged.', $source);
                return 1;
            } else {
                $self->fatal_message('Cannot stage results. Multiple differing outputs with identical names: %s', $file);
            }
        }

        Genome::Sys->move($source, $destination);
    }

    return 1;
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
    my $results_dir = shift;

    Genome::Sys->remove_directory_tree($tmp_dir);

    if ($results_dir) {
        for my $dir (glob("$results_dir/tmp*")) {
            if (-d $dir) {
                Genome::Sys->remove_directory_tree($dir);
            }
        }
    }

    return 1;
}

1;
