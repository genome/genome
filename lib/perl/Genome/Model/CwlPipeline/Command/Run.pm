package Genome::Model::CwlPipeline::Command::Run;

use strict;
use warnings;

use Cwd;
use Data::Dumper;
use File::Spec;
use Genome;
use Genome::Utility::File::Mode qw();
use Genome::Utility::Text qw();
use YAML;
use JSON qw(to_json from_json);
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
    } elsif ($cwl_runner eq 'cromwell_gcp') {
        $self->run_cromwell_gcp($yaml, $tmp_dir, $results_dir);
        $self->cleanup($tmp_dir);
    } elsif ($cwl_runner =~ m'^toil') {
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

    my $cwl_runner = Genome::Config::get('cwl_runner');
    my($toil, $version) = split(" ", $cwl_runner);
    if ($toil ne 'toil') {
        $self->fatal_message('Called run_toil with non-toil runner: %s', $toil);
    }
    my $wrapper = Genome::Model::CwlPipeline::Runner::Toil->cwl_runner_wrapper_for_version($version);

    my $build = $self->build;
    my $model = $build->model;
    my $jobstore_dir = File::Spec->join($tmp_dir, 'jobstore');
    my $work_dir = File::Spec->join($tmp_dir, 'work');
    my $toil_tmp_output_dir = File::Spec->join($tmp_dir, 'toil-tmp');
    my $log_file = File::Spec->join($build->log_directory, 'toil.log');

    my @restart;
    if (-e $jobstore_dir) {
        push @restart, '--restart';
    }
    Genome::Sys->create_directory($work_dir);

    my $primary_docker_image = $model->primary_docker_image;
    local $ENV{LSB_SUB_ADDITIONAL} = $primary_docker_image;

    my $default_queue = Genome::Config::get('lsf_queue_build_worker_alt');
    local $ENV{LSB_DEFAULTQUEUE} = $default_queue;

    my $lsf_group = Genome::Config::get('lsf_user_group');
    local $ENV{LSB_SUB_USER_GROUP} = $lsf_group;

    local $ENV{LSF_DOCKER_PRESERVE_ENVIRONMENT} = 'false';

    local $ENV{TOIL_CHECK_ENV} = 'True';

    #toil relies on reading the output from bsub
    delete local $ENV{BSUB_QUIET};

    Genome::Sys->shellcmd(
        cmd => [
            '/bin/bash',
            $wrapper,
            '--disableCaching', '--logLevel=DEBUG',
            @restart,
            '--workDir', $work_dir,
            '--tmp-outdir-prefix', $toil_tmp_output_dir,
            '--jobStore', $jobstore_dir,
            "--logFile=$log_file",
            "--outdir=$results_dir",
            "--no-container", #with patched toil, this allows per-job containers
            '--batchSystem', 'lsf',
            '--writeLogs', $build->log_directory,
            '--writeLogsFromAllJobs',
            '--maxLogFileSize', 1_000_000_000, #bytes
            '--bypass-file-store',
            '--stats',
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

    my @jar_cmdline = Genome::Cromwell->cromwell_jar_cmdline($config_file);

    Genome::Sys->shellcmd(
        cmd => [
            @jar_cmdline,
            'run',
            '-t', $wf_type,
            '-l', $labels_file,
            '-i', $yaml,
            $main_workflow_file,
        ],
        input_files => [$main_workflow_file, $yaml],
    );

    my $guard = Genome::Cromwell->spawn_local_server($config_file);
    $self->_stage_cromwell_outputs($results_dir);
}

sub run_cromwell_gcp {
    my $self = shift;
    my $yaml = shift;
    my $tmp_dir = shift;
    my $results_dir = shift;

    my $logdir = $self->build->log_directory;
    my $data_dir = $self->build->data_directory;
    my $main_workflow_file = $self->build->model->main_workflow_file;
    my $build_id = $self->build->id;

    my $bucket = Genome::Config::get('cromwell_gcp_bucket');
    my $cromwell_gcp_project = Genome::Config::get('cromwell_gcp_project');
    my $cromwell_server_memory_gb = Genome::Config::get('cromwell_gcp_server_memory_gb');
    my $cromwell_service_account = Genome::Config::get('cromwell_gcp_service_account');
    my $cromwell_subnet = Genome::Config::get('cromwell_gcp_subnet');
    my $queue = Genome::Config::get('lsf_queue_build_worker');
    my $user_group = Genome::Config::get('lsf_user_group');

    my $poll_interval_seconds = 300;

    #
    # Cloudize workflow
    #
    my $cloud_yaml = $yaml; $cloud_yaml =~ s/.ya?ml/_cloud.json/;
    my $lsb_sub_guard = Genome::Config::set_env('lsb_sub_additional', 'docker(mgibio/cloudize-workflow:1.2.2)');
    delete local $ENV{BOTO_CONFIG};

    Genome::Sys::LSF::bsub::bsub(
        queue => $queue,
        user_group => $user_group,
        resource_string => 'rusage[mem=512M:internet2_upload_mbps=500]',
        wait_for_completion => 1,
        log_file => File::Spec->join($logdir, '00_cloudize_workflow.log'),
        cmd => [
            "python3",
            "/opt/scripts/cloudize-workflow.py",
            $bucket,
            $main_workflow_file,
            $yaml,
            "--output=$cloud_yaml"] );

    # Zip dependencies
    my $deps_zip_path = File::Spec->join($data_dir, 'deps.zip');
    my $deps_zip_url = "gs://$bucket/build.$build_id/deps.zip";
    my $prev_dir = getcwd;
    my(undef, $deps_dir, undef) = File::Spec->splitpath($main_workflow_file);

    chdir($deps_dir);
    Genome::Sys->shellcmd(
        cmd => ['zip', '-r', $deps_zip_path, '.'],
        redirect_stdout => '/dev/null'
        );
    chdir($prev_dir);

    # Upload zip file
    Genome::Sys::LSF::bsub::bsub(
        queue => $queue,
        user_group => $user_group,
        resource_string => 'rusage[internet2_download_mbps=500]',
        wait_for_completion => 1,
        log_file => File::Spec->join($logdir, '01_upload_zip.log'),
        cmd => ['gsutil', 'cp', '-n', $deps_zip_path, $deps_zip_url]
        );

    #
    # Generate files for run
    #
    my $conf_file = $self->_generate_cromwell_config_gcp($tmp_dir);
    my $options_file = $self->_generate_workflow_options_gcp;
    my $labels_file = $self->_generate_cromwell_labels;

    #
    # Do the run
    #
    $self->status_message("Files generated. Starting VM.");
    Genome::Sys::LSF::bsub::bsub(
        queue => $queue,
        user_group => $user_group,
        wait_for_completion => 1,
        log_file => File::Spec->join($logdir, '02_vm_start.log'),
        cmd => [
            "sh", "/opt/gms/start.sh",
            "--build", $build_id,
            "--cromwell-conf", $conf_file,
            "--service-account", $cromwell_service_account,
            "--workflow-definition", $main_workflow_file,
            "--workflow-inputs", $cloud_yaml,
            "--workflow-options", $options_file,
            "--deps-zip", $deps_zip_url,
            "--bucket", $bucket,
            "--subnet", $cromwell_subnet,
            "--project", $cromwell_gcp_project,
            "--memory-gb", $cromwell_server_memory_gb,
            "--tmp-dir", $tmp_dir
        ] );

    my $result;
    # Wait for instance VM to terminate itself
    $self->status_message("VM started. Polling every $poll_interval_seconds seconds.");
    do {
        sleep $poll_interval_seconds;
        $result = system("gcloud compute instances describe build-$build_id --zone us-central1-c > /dev/null");
        $self->status_message("Polled VM and got result $result");
    } while ($result == 0);
    $self->status_message("Polling done. Pulling artifacts.\n");

    # Pull build directory
    Genome::Sys::LSF::bsub::bsub(
        queue => $queue,
        user_group => $user_group,
        resource_string => 'rusage[internet2_download_mbps=500]',
        wait_for_completion => 1,
        log_file => File::Spec->join($logdir, '03_pull_dir.log'),
        cmd => ['gsutil', 'cp', '-r', '-n', "gs://$bucket/build.$build_id/*", $data_dir]
        );
    $self->status_message("Pulled artifacts. Pulling outputs.");

    # Fetch outputs files
    my $outputs_json = File::Spec->join($data_dir, 'outputs.json');
    # Pull output files
    if (-e $outputs_json) {
        Genome::Sys::LSF::bsub::bsub(
            queue => $queue,
            user_group => $user_group,
            resource_string => 'rusage[mem=512M:internet2_download_mbps=500]',
            wait_for_completion => 1,
            log_file => File::Spec->join($logdir, '04_pull_outputs.log'),
            cmd => [
                'python3', '/opt/scripts/pull_outputs.py',
                "--outputs-dir=$results_dir",
                "--outputs-file=$outputs_json"
            ] );
    } else {
        $self->fatal_message("Build did not generate output files. See $logdir");
    }
    $self->status_message("Finished Google Cloud run of workflow.\n" .
                          "Results in $results_dir \n" .
                          "Compute instance logs and workflow timing in $data_dir \n" .
                          "Logs for bsubs at $logdir" );
}

sub _fetch_cromwell_log {
    my $self = shift;
    my $workflow_id = shift;
    my $workflow_options = shift;
    my $logdir = shift;

    my $workflow_opts = from_json(Genome::Sys->read_file($workflow_options));
    my $final_workflow_log_dir = $workflow_opts->{final_workflow_log_dir};

    Genome::Sys->shellcmd(
        cmd => [
            '/usr/bin/python3',
            '/usr/bin/gsutil/gsutil', 'cp',
            File::Spec->join($final_workflow_log_dir, "workflow." . $workflow_id . ".log"),
            File::Spec->join($logdir, $workflow_id . ".log")
        ]);
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

sub _generate_workflow_options_gcp {
    my $self = shift;

    my $build_id = $self->build->id;
    my $data_dir = $self->build->data_directory;
    my $cromwell_gcp_bucket = Genome::Config::get('cromwell_gcp_bucket');
    my $options_file = File::Spec->join($data_dir, 'gcp_workflow_options.json');

    return $options_file if -e $options_file;

    my $contents = <<"EOCONFIG"
{
    "final_workflow_log_dir": "gs://$cromwell_gcp_bucket/build.$build_id/logs",
    "final_call_logs_dir": "gs://$cromwell_gcp_bucket/build.$build_id/logs",
    "use_relative_output_paths": false,
    "default_runtime_attributes": {
        "preemtible": 1
    }
}
EOCONFIG
        ;
    Genome::Sys->write_file($options_file, $contents);
    return $options_file;
}

sub _generate_cromwell_config_gcp {
    my $self = shift;
    my $tmp_dir = shift;

    my $data_dir = $self->build->data_directory;
    my $config_file = File::Spec->join($data_dir, 'cromwell_gcp.config');
    return $config_file if -e $config_file;

    my $cromwell_gcp_project = Genome::Config::get('cromwell_gcp_project');
    my $cromwell_gcp_service_account = Genome::Config::get('cromwell_gcp_service_account');
    my $cromwell_gcp_bucket = Genome::Config::get('cromwell_gcp_bucket');

    my $build_id = $self->build->id;

    my $config = <<'EOCONFIG'
include required(classpath("application"))

google {
  application-name = "cromwell"
  auths = [
    {
      name = "application-default"
      scheme = "application_default"
    }
  ]
}

engine.filesystems {
  gcs.auth = "application-default"
  local.enabled = true
}

backend.default = "default"
backend.providers.default {
  actor-factory = "cromwell.backend.google.pipelines.v2beta.PipelinesApiLifecycleActorFactory"
  config {
    genomics {
      auth = "application-default"
      endpoint-url = "https://lifesciences.googleapis.com/"
      location = "us-central1"
    }
    filesystems {
      gcs {
        auth = "application-default"
        caching.duplication-strategy = "reference"
      }
    }
    include "papi_v2_reference_image_manifest.conf"
  }
}
EOCONFIG
        ;

    $config .= <<"EOCONFIG"
backend.providers.default.config {
  project = $cromwell_gcp_project
  root = "gs://$cromwell_gcp_bucket/build.$build_id"
  genomics.compute-service-account = "$cromwell_gcp_service_account"
  filesystems.gcs.project = $cromwell_gcp_project
}
EOCONFIG
        ;

    Genome::Sys->write_file($config_file, $config);
    return $config_file;
}

sub _generate_cromwell_config {
    my $self = shift;
    my $tmp_dir = shift;
    my $build = $self->build;

    my $data_dir = $build->data_directory;
    my $config_file = File::Spec->join($data_dir,'cromwell.config');
    return $config_file if -e $config_file;

    my $log_dir = $build->log_directory;

    my $primary_docker_image = $build->model->primary_docker_image;

    my $default_queue = Genome::Config::get('lsf_queue_build_worker_alt');

    my $job_group = join('/',
        Genome::Config::get('lsf_job_group'),
        Genome::Sys->username,
        'workers'
    );
    my $user_group = Genome::Config::get('lsf_user_group');

    my $server = Genome::Config::get('cromwell_server');
    my $auth = Genome::Config::get('cromwell_auth');
    my $user = Genome::Config::get('cromwell_user');

    my $docker_volumes = Genome::Config::get('docker_volumes');
    if ($ENV{LSF_DOCKER_VOLUMES}) {
        if ($docker_volumes) {
            $docker_volumes = $ENV{LSF_DOCKER_VOLUMES} . ' ' . $docker_volumes;
        } else {
            $docker_volumes = $ENV{LSF_DOCKER_VOLUMES};
        }
    }

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
        -G '$user_group' \\
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
        -G '$user_group' \\
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
        kill-docker = "bkill ${job_id}"
        check-alive = "bjobs -noheader -o stat ${job_id} | /bin/grep 'PEND\\|RUN'"
        job-id-regex = "Job <(\\d+)>.*"
EOCONFIG
;
    $config .= <<EOCONFIG
        root = "$tmp_dir/cromwell-executions"
EOCONFIG
;
    if(Genome::Config::get('cromwell_call_caching')) {
        $config .= <<'EOCONFIG'
        exit-code-timeout-seconds = 600

        filesysytems {
          local {
            caching {
              duplication-strategy: [
                "hard-link", "soft-link", "copy"
              ]
              hashing-strategy: "xxh64"
              fingerprint-size: 10485760
              check-sibling-md5: false
            }
          }
        }
EOCONFIG
;
    }

    $config .= <<EOCONFIG
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

        #shell out so this is saved immediately for the benefit of `genome model build view` while this build runs.
        Genome::Sys->shellcmd(
            cmd => [qw(genome model build add-note --header-text=hsqldb_server_file), "--body-text=$dbfile_location", $build->id],
        );

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

    if (Genome::Config::get('cromwell_call_caching')) {
        $config .= <<EOCONFIG
call-caching {
  enabled = true
  invalidate-bad-cache-results = true
}
EOCONFIG
;
    }

    Genome::Sys->write_file($config_file, $config);
    return $config_file;
}

sub _generate_cromwell_labels {
    my $self = shift;
    my $build = $self->build;

    my $labels_file = File::Spec->join($build->data_directory, 'cromwell.labels');

    return $labels_file if -e $labels_file;

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

    my $results = Genome::Cromwell->query( [{ label => 'build:' . $build->id, status => 'Succeeded' }] );
    if ($results->{totalResultsCount} < 1) {
        $build->fatal_message('Failed to find workflow.  Got: %s', scalar(Data::Dumper::Dumper($results)));
    }

    my $workflows = $results->{results};
    if (@$workflows > 1) {
        $workflows = [sort { $b->{end} cmp $a->{end} } @$workflows];
    }

    my $workflow_id = $workflows->[0]->{id};
    my $output_result = Genome::Cromwell->outputs($workflow_id);

    my $outputs = $output_result->{outputs};

    my $prefix = $self->_determine_output_prefix;
    for my $output_name (keys %$outputs) {
        my $info = $outputs->{$output_name};

        $self->_stage_cromwell_output($results_dir, $info, $prefix);
    }

    $self->_generate_timing_report($workflow_id);

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

sub _generate_timing_report {
    my $self = shift;
    my $workflow_id = shift;

    my $data_dir = $self->build->data_directory;
    my $timing_file = File::Spec->join($data_dir, 'timing.html');

    my $report = Genome::Cromwell->timing($workflow_id);
    Genome::Sys->write_file($timing_file, $report);

    return $timing_file;
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
