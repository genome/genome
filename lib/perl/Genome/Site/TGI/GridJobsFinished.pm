package Genome::Site::TGI::GridJobsFinished;

use strict;
use warnings;
use Genome;
use POSIX;

class Genome::Site::TGI::GridJobsFinished {
    table_name => "GRID_JOBS_FINISHED",
    id_by => [
        bjob_id => {
            is => 'Number',
        },
    ],
    has => [
        job_id => {
            is => 'Number',
            column_name => 'jobid',
        },
        index_id => {
            is => 'Number',
            column_name => 'indexid',
        },
        cluster_id => {
            is => 'Number',
            column_name => 'clusterid',
        },
        options => {
            is => 'Number',
        },
        options2 => {
            is => 'Number',
        },
        options3 => {
            is => 'Number',
        },
        stat_changes => {
            is => 'Number',
        },
        flapping_logged => {
            is => 'Number',
        },
        exit_status => {
            is => 'Number',
            column_name => 'exitstatus',
        },
        exec_uid => {
            is => 'Number',
            column_name => 'execuid',
        },
        job_priority => {
            is => 'Number',
            column_name => 'jobpriority',
        },
        job_pid => {
            is => 'Number',
            column_name => 'jobpid',
        },
        user_priority => {
            is => 'Number',
            column_name => 'userpriority',
        },
        cpu_used => {
            is => 'Number',
        },
        u_time => {
            is => 'Number',
            column_name => 'utime',
        },
        s_time => {
            is => 'Number',
            column_name => 'stime',
        },
        efficiency => {
            is => 'Number',
        },
        efficiency_logged => {
            is => 'Number',
            column_name => 'effic_logged',
        },
        num_pids => {
            is => 'Number',
            column_name => 'numpids',
        },
        num_pgids => {
            is => 'Number',
            column_name => 'numpgids',
        },
        num_threads => {
            is => 'Number',
            column_name => 'numthreads',
        },
        pid_alarm_logged => {
            is => 'Number',
        },
        num_nodes => {
            is => 'Number',
        },
        num_cpus => {
            is => 'Number',
        },
        max_num_processors => {
            is => 'Number',
            column_name => 'maxnumprocessors',
        },
        submit_time => {
            is => 'DateTime',
        },
        reserve_time => {
            is => 'DateTime',
            column_name => 'reservetime',
        },
        predicted_start_time => {
            is => 'DateTime',
            column_name => 'predictedstarttime',
        },
        start_time => {
            is => 'DateTime',
        },
        end_time => {
            is => 'DateTime',
        },
        begin_time => {
            is => 'DateTime',
            column_name => 'begintime',
        },
        term_time => {
            is => 'DateTime',
            column_name => 'termtime',
        },
        completion_time => {
            is => 'Number',
        },
        pend_time => {
            is => 'Number',
        },
        psusp_time => {
            is => 'Number',
        },
        run_time => {
            is => 'Number',
        },
        ususp_time => {
            is => 'Number',
        },
        ssusp_time => {
            is => 'Number',
        },
        unknown_time => {
            is => 'Number',
            column_name => 'unkwn_time',
        },
        rlimit_max_cpu => {
            is => 'Number',
        },
        rlimit_max_wallt => {
            is => 'Number',
        },
        rlimit_max_swap => {
            is => 'Number',
        },
        rlimit_max_fsize => {
            is => 'Number',
        },
        rlimit_max_data => {
            is => 'Number',
        },
        rlimit_max_stack => {
            is => 'Number',
        },
        rlimit_max_core => {
            is => 'Number',
        },
        rlimit_max_rss => {
            is => 'Number',
        },
        job_start_logged => {
            is => 'Number',
        },
        job_end_logged => {
            is => 'Number',
        },
        job_scan_logged => {
            is => 'Number',
        },
        last_updated => {
            is => 'DateTime',
        },
    ],
    has_optional => [
        bjob_user => {
            is => 'Text',
        },
        'stat' => {
            is => 'Text',
        },
        prev_stat => {
            is => 'Text',
        },
        pend_reasons => {
            is => 'Text',
            column_name => 'pendreasons',
        },
        queue => {
            is => 'Text',
        },
        nice => {
            is => 'Text',
        },
        from_host => {
            is => 'Text',
        },
        exec_host => {
            is => 'Text',
        },
        login_shell => {
            is => 'Text',
            column_name => 'loginshell',
        },
        exec_home => {
            is => 'Text',
            column_name => 'exechome',
        },
        exec_cwd => {
            is => 'Text',
            column_name => 'execcwd',
        },
        cwd => {
            is => 'Text',
        },
        post_exec_cmd => {
            is => 'Text',
            column_name => 'postexeccmd',
        },
        exec_user_name => {
            is => 'Text',
            column_name => 'execusername',
        },
        mail_user => {
            is => 'Text',
            column_name => 'mailuser',
        },
        job_name => {
            is => 'Text',
            column_name => 'jobname',
        },
        project_name => {
            is => 'Text',
            column_name => 'projectname',
        },
        parent_group => {
            is => 'Text',
            column_name => 'parentgroup',
        },
        sla => {
            is => 'Text',
        },
        job_group => {
            is => 'Text',
            column_name => 'jobgroup',
        },
        license_project => {
            is => 'Text',
            column_name => 'licenseproject',
        },
        command => {
            is => 'Text',
        },
        new_command => {
            is => 'Text',
            column_name => 'newcommand',
        },
        in_file => {
            is => 'Text',
            column_name => 'infile',
        },
        out_file => {
            is => 'Text',
            column_name => 'outfile',
        },
        err_file => {
            is => 'Text',
            column_name => 'errfile',
        },
        pre_exec_cmd => {
            is => 'Text',
            column_name => 'preexeccmd',
        },
        resource_requirements => {
            is => 'Text',
            column_name => 'res_requirements',
        },
        depend_condition => {
            is => 'Text',
            column_name => 'dependcond',
        },
        memory_used => {
            is => 'Text',
            column_name => 'mem_used',
            default => 0,
        },
        max_memory => {
            is => 'Number',
            default => 0,
        },
        max_swap => {
            is => 'Number',
        },
        memory_requested => {
            is => 'Number',
            column_name => 'mem_requested',
            default => 0,
        },
        memory_requested_operation => {
            is => 'Text',
            column_name => 'mem_requested_oper',
        },
        host_spec => {
            is => 'Text',
            column_name => 'hostspec',
        },
        user_group => {
            is => 'Text',
            column_name => 'usergroup',
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub cpu_instance_size {
    return 1;
}

sub memory_instance_size {
    return 4194000; # 4 GB
}

sub cpu_instances_used {
    my $self = shift;
    return POSIX::ceil($self->num_cpus / $self->cpu_instance_size);
}

# Memory requested can be 0 if the user didn't actually request any, or they
# could have under-requested. In either case, the actual memory usage is stored.
sub actual_memory_used {
    my $self = shift;
    return List::Util::max($self->memory_requested, $self->max_memory);
}

sub memory_instances_used {
    my $self = shift;
    return POSIX::ceil($self->actual_memory_used / $self->memory_instance_size);
}

sub total_instances_used {
    my $self = shift;
    return List::Util::max($self->cpu_instances_used, $self->memory_instances_used);
}

sub wallclock_time_seconds {
    my $self = shift;
    return $self->run_time;
}

sub wallclock_time_hours {
    my $self = shift;
    return $self->wallclock_time_seconds / (60 * 60);
}

sub total_instance_hours_used {
    my $self = shift;
    return $self->total_instances_used * $self->wallclock_time_hours;
}

1;
