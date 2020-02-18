package Genome::Model::CwlPipeline::Command::Run;

use strict;
use warnings;

use Data::Dumper;
use File::Spec;
use Genome;
use Genome::Utility::File::Mode qw();
use Genome::Utility::Text qw();
use YAML;
use JSON qw(to_json);
use File::Compare qw();
use File::Copy::Recursive qw();

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
    doc => 'wrapper command to run the workflow for a build'
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
    $self->process_outputs($results_dir);

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

    my $tmp_dir = $build->get_or_create_scratch_directory;

    my $data_directory = $build->data_directory;
    my $results_dir = File::Spec->join($data_directory, 'results');
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

    my $main_workflow_file = $model->main_workflow_file;

    my $wf_type = $self->_determine_workflow_type($main_workflow_file);

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
            '-t', $wf_type,
            '-l', $labels_file,
            '-i', $yaml,
            $main_workflow_file,
        ],
        input_files => [$main_workflow_file, $yaml],
    );

    $self->_stage_cromwell_outputs($results_dir);
}

sub _determine_workflow_type {
    my $self = shift;
    my $main_workflow_file = shift;

    if($main_workflow_file =~ /\.wdl$/) {
        return 'wdl';
    } else {
        return 'cwl';
    }
}

sub _generate_cromwell_config {
    my $self = shift;
    my $tmp_dir = shift;

    my $build = $self->build;
    my $log_dir = $build->log_directory;

    my $primary_docker_image = $build->model->primary_docker_image;

    my $default_queue = Genome::Config::get('lsf_queue_build_worker_alt');

    my $job_group = join('/',
        Genome::Config::get('lsf_job_group'),
        Genome::Sys->username,
        'workers'
    );

    my $server = Genome::Config::get('cromwell_server');
    my $auth = Genome::Config::get('cromwell_auth');
    my $user = Genome::Config::get('cromwell_user');

    my $docker_volumes = Genome::Config::get('docker_volumes');

    my $data_dir = $build->data_directory;
    my $config_file = File::Spec->join($data_dir,'cromwell.config');

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
        Int memory_kb = 4096000
        Int memory_mb = 4096
        String? docker
        """

        submit = """
EOCONFIG
    ;
    if ($docker_volumes) {
        $config .= <<EOCONFIG
        LSF_DOCKER_VOLUMES='$docker_volumes' \
EOCONFIG
        ;
    }
    $config .= <<'EOCONFIG'
        LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
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
        -g '$job_group' \\
EOCONFIG
        ;
    $config .= <<'EOCONFIG'
        -M ${memory_kb} \
        -n ${cpu} \
        -R "span[hosts=1] select[mem>${memory_mb}] rusage[mem=${memory_mb}]" \
        /bin/bash ${script}
        """

        submit-docker = """
EOCONFIG
    ;

    my $vol = '${cwd}:${docker_cwd}';
    if ($docker_volumes) {
        $vol = join(' ', $vol, $docker_volumes);
    }
    $config .= <<EOCONFIG
        LSF_DOCKER_VOLUMES='$vol' \
EOCONFIG
    ;
    $config .= <<'EOCONFIG'
        LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
        bsub \
        -J ${job_name} \
        -cwd ${cwd} \
EOCONFIG
    ;
    $config .= <<EOCONFIG
        -o /dev/null \\
        -e $log_dir/cromwell-%J.err \\
        -q '$default_queue' \\
        -g '$job_group' \\
EOCONFIG
        ;
    $config .= <<'EOCONFIG'
        -a "docker(${docker})" \
        -M ${memory_kb} \
        -n ${cpu} \
        -R "span[hosts=1] select[mem>${memory_mb}] rusage[mem=${memory_mb}]" \
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
EOCONFIG
;
    if ($server =~ /^mysql:/) {
        $config .= <<EOCONFIG
database {
  profile = "slick.jdbc.MySQLProfile\$"
  db {
    driver = "com.mysql.jdbc.Driver"
    url = "jdbc:$server"
    user = "$user"
    password = "$auth"
    connectionTimeout = 30000
    numThreads = 5
  }
}
EOCONFIG
;
    } elsif ($server =~ /^hsqldb:/) {
        my $dbfile_location;
        if ($server =~ /;/) {
            $self->fatal_message('Cannot currently handle hsqldb server string with semicolons. Got: %s', $server);
        } elsif ($server eq 'hsqldb:tmp') {
            $dbfile_location = "$tmp_dir/cromwell-db/cromwell-db";
            $server = "hsqldb:file:$dbfile_location";
            $self->debug_message('Using temporary hsqldb location: %s', $server);
        } elsif ($server eq 'hsqldb:build') {
            $dbfile_location = "$data_dir/cromwell-db/cromwell-db";
            $server = "hsqldb:file:$dbfile_location";
            $self->debug_message('Using build hsqldb location: %s', $server);
        } else {
            ($dbfile_location) = $server =~ /hsqldb:file:([^:]+)/;
            unless ($dbfile_location) {
                $self->fatal_message('Could not parse hsqldb file location from server string. Expected "hsqldb:tmp", "hsqldb:build", or "hsqldb:file:/path/to/cromwell-db". Got: %s', $server);
            }
            $self->debug_message('Using supplied hsqldb location: %s', $dbfile_location);
        }

        my $note = $build->add_note(header_text => 'hsqldb_server_file', body_text => $dbfile_location);
        $note->body_text($dbfile_location); #ignore any system-generated sudo message for this note.

        $config .= <<EOCONFIG
database {
  profile = "slick.jdbc.HsqldbProfile\$"
  db {
    driver = "org.hsqldb.jdbcDriver"
    url = """
    jdbc:$server;
    shutdown=false;
    hsqldb.default_table_type=cached;hsqldb.tx=mvcc;
    hsqldb.result_max_memory_rows=10000;
    hsqldb.large_data=true;
    hsqldb.applog=1;
    hsqldb.lob_compressed=true;
    hsqldb.script_format=3
    """
    connectionTimeout = 120000
    numThreads = 1
   }
}
EOCONFIG
;
    } else {
        $self->fatal_message('Expected mysql or hsqldb cromwell server url but got: %s', $server);
    }

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
        $self->fatal_message('Failed to find workflow.  Got: %s', scalar(Data::Dumper::dumper($results)));
    }

    my $workflow_id = $results->{results}->[0]->{id};
    my $output_result = Genome::Cromwell->outputs($workflow_id);

    my $outputs = $output_result->{outputs};

    my $prefix = $self->_determine_output_prefix;
    for my $output_name (keys %$outputs) {
        my $info = $outputs->{$output_name};

        $self->_stage_cromwell_output($results_dir, $info, $prefix);
    }

    return 1;
}

sub _determine_output_prefix {
    my $self = shift;

    my $build = $self->build;

    my $prefix;
    my (@prefix_input) = grep { $_->name eq 'output_prefix' } $build->inputs;
    if (@prefix_input) {
        $prefix = join('_', sort map { $_->value_id } @prefix_input);
        $prefix = Genome::Utility::Text::sanitize_string_for_filesystem($prefix);
    }

    return $prefix;
}

sub _stage_cromwell_output {
    my $self = shift;
    my $results_dir = shift;
    my $item = shift;
    my $prefix = shift;

    if (ref $item eq 'ARRAY') {
        map $self->_stage_cromwell_output($results_dir, $_, $prefix), @$item;
        return 1;
    }

    return unless ref $item eq 'HASH';
    my $location = $item->{location};
    return unless $location;

    my @secondary = map { $_->{location} } @{ $item->{secondaryFiles} };

    for my $source ($location, @secondary) {
        my (undef, $dir, $file) = File::Spec->splitpath($source);

        if ($prefix) {
            $file = join('-', $prefix, $file);
        }

        my $destination = File::Spec->join($results_dir, $file);
        if (-e $destination) {
            if( File::Compare::compare($source, $destination) == 0 ) {
                $self->warning_message('Skipping staging of %s--an identical file was already staged.', $source);
                return 1;
            } else {
                $self->fatal_message('Cannot stage results. Multiple differing outputs with identical names: %s', $file);
            }
        }

        if (-l $source) {
            if ($prefix) {
                my $target = readlink $source;
                symlink(join('-', $prefix, $target), $destination);
            } else {
                File::Copy::Recursive::fmove($source, $destination);
            }
        } elsif (-d $source) {
            File::Copy::Recursive::dirmove($source, $destination);
        } else {
            Genome::Sys->move($source, $destination);
        }
    }

    return 1;
}

sub preserve_results {
    my $self = shift;
    my $results_dir = shift;

    for my $file ($results_dir, glob("$results_dir/*")) {
        next unless -e $file; #sometimes workflows produce dangling symlinks, yay!
        Genome::Utility::File::Mode::mode($file)->rm_all_writable;
    }

    return 1;
}

sub cleanup {
    my $self = shift;
    my $tmp_dir = shift;
    my $results_dir = shift;

    my $build = $self->build;
    $build->cleanup_scratch_directory;

    if ($results_dir) {
        for my $dir (glob("$results_dir/tmp*")) {
            if (-d $dir) {
                Genome::Sys->remove_directory_tree($dir);
            }
        }
    }

    return 1;
}

sub process_outputs {
    my $self = shift;
    my $results_dir = shift;

    my $build = $self->build;
    my $pp = $build->processing_profile;
    my $da = Genome::Disk::Allocation->get(owner => $pp);
    my $dir = $da->absolute_path;

    my $postprocessing_script = File::Spec->join($dir, 'process_outputs.pl');
    if (-e $postprocessing_script) {
        Genome::Sys->shellcmd(
            cmd => [$^X, $postprocessing_script, $build->id, $results_dir],
        );
    }

    return 1;
}

1;
