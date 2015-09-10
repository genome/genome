package Genome::Utility::LSFResourceParser;

use strict;
use warnings FATAL => qw(all);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(parse_lsf_params construct_lsf_param_string);

use Getopt::Long qw(GetOptionsFromString);
use IO::String;

sub parse_lsf_params {
    my ($lsf_param_string) = @_;

    my ($parse_ok, $lsf_params, $message) = _parse_lsf_params($lsf_param_string);

    unless ($parse_ok) {
        ($parse_ok, $lsf_params) = _parse_lsf_params(qq{-R '$lsf_param_string'});
        unless ($parse_ok) {
            die sprintf("Failed to parse lsf options: %s\n%s\n",
                $lsf_param_string, $message);
        }
    }

    return $lsf_params;
}

sub _parse_lsf_params {
    my ($lsf_param_string) = @_;

    my %lsf_params = ('options' => {}, 'rLimits' => {});
    my ($parse_ok, $message) = _get_options_from_string($lsf_param_string,
        _create_getopt_specs($lsf_params{options}, _valid_options()),
        _create_getopt_specs($lsf_params{rLimits}, _valid_rlimits()));

    return ($parse_ok, \%lsf_params, $message);
}

sub _get_options_from_string {
    my $getopt_fh = IO::String->new;
    local *STDERR = $getopt_fh;
    my $ret = GetOptionsFromString(@_);
    return ($ret, ${$getopt_fh->string_ref});
}

sub _create_getopt_specs {
    my ($destination, %option) = @_;

    my @getopt_specs;
    while (my ($bsub_option, $submit_field) = each %option) {
        push @getopt_specs, "$bsub_option=s" => sub {
            my ($option_name, $option_value) = @_;
            $destination->{$submit_field} = $option_value;
        };
    }

    return @getopt_specs;
}

sub _create_option_spec {
    my ($option, $dictionary) = @_;
    return ("$option=s" => sub { my ($option_name, $value, )});
}

my %valid_options = (
    'b' => 'beginTime',
    'e' => 'errFile',
    'g' => 'group',
    'i' => 'inFile',
    'J' => 'jobName',
    'u' => 'mail_user',
    'n' => 'numProcessors',  # ,maxNumProcessors
    'o' => 'outFile',
    'E' => 'preExecCmd',
    'Ep' => 'postExecCmd',
    'P' => 'projectName',
    'q' => 'queue',
    'R' => 'resReq',
    't' => 'termTime',
);
my %option_lookup = _reverse_hash(%valid_options);
sub _valid_options { %valid_options }
sub _option_lookup { %option_lookup }


my %valid_rlimits = (
    'c' => 'cpuTime',
    'M' => 'RSS',
    'F' => 'openFiles',
    'S' => 'stack',
);
my %rlimit_lookup = _reverse_hash(%valid_rlimits);
sub _valid_rlimits { %valid_rlimits }
sub _rlimit_lookup { %rlimit_lookup }



sub _reverse_hash {
    my %hash = @_;
    my %reversed_hash;
    while (my ($k,$v) = each %hash) {
        $reversed_hash{$v} = $k;
    }
    return %reversed_hash;
}

sub construct_lsf_param_string {
    my ($lsf_params) = @_;

    my $lsf_param_string = join(q{ },
        _construct_lsf_param_string(
            $lsf_params->{options}, _option_lookup()),
        _construct_lsf_param_string(
            $lsf_params->{rLimits}, _rlimit_lookup()));

    return $lsf_param_string;
}

sub _construct_lsf_param_string {
    my ($source, %param_lookup) = @_;

    my @params;
    while (my ($submit_field, $value) = each %$source) {
        my $param = _lookup_param($submit_field, \%param_lookup);
        unless (defined $param) {
            die sprintf(
                'Failed to construct lsf params.  Unkown submit field: %s',
                $submit_field);
        }
        push @params, "-$param", "'$value'";
    }

    if (@params) {
        return join(q{ }, @params);
    }
    else {
        return;
    }
}

sub _lookup_param {
    my ($submit_field, $lookup) = @_;
    return exists($lookup->{$submit_field}) && $lookup->{$submit_field};
}


1;
__END__


This table is copied from "Platform LSF Programmer's Guide" for LSF 8.0 pages 67-68
bsub Option                       submitField                                 options
-------------------------------   -------------------------------             ---
-J job_name_spec                  jobName                                     SUB_JOB_NAME
-q queue_name                     queue                                       SUB_QUEUE
-m host_name[+[pref_level]]       askedHosts                                  SUB_HOST
-n min_proc[,max_proc]            numProcessors, maxNumProcessors
-R res_req                        resReq                                      SUB_RES_REQ
-c cpu_limit[/host_spec]          rlimits[LSF_RLIMIT_CPU] / hostSpec **       SUB_HOST_SPEC (if host_spec is specified)
-W run_limit[/host_spec]          rlimits[LSF_RLIMIT_RUN] / hostSpec**        SUB_HOST_SPEC (if host_spec is specified)
-F file_limit                     rlimits[LSF_RLIMIT_FSIZE]**
-M mem_limit                      rlimits[LSF_RLIMIT_RSS]**
-D data_limit                     rlimits[LSF_RLIMIT_DATA]**
-S stack_limit                    rlimits[LSF_RLIMIT_STACK**
-C core_limit                     rlimits[LSF_RLIMIT_CORE]**
-k "chkpnt_dir [chkpnt_period]"   chkpntDir, chkpntPeriod                     SUB_CHKPNT_DIR, SUB_CHKPNT_DIR (if chkpntPeriod is specified)
-w depend_cond                    dependCond                                  SUB_DEPEND_COND
-b begin_time                     beginTime
-t term_time                      TermTime
-i in_file                        inFile                                      SUB_IN_FILE
-o out_file                       outFile                                     SUB_OUT_FILE
-e err_file                       errFile                                     SUB_ERR_FILE
-u mail_user                      mailUser                                    SUB_MAIL_USER
-f "lfile op [rfile]"             xf
-E "pre_exec_cmd [arg]"           preExecCmd                                  SUB_PRE_EXEC
-L login_shell                    loginShell                                  SUB_LOGIN_SHELL
-P project_name                   projectName                                 SUB_PROJECT_NAME
-G user_group                     userGroup                                   SUB_USER_GROUP
-H                                                                            SUB2_HOLD*
-x                                                                            SUB_EXCLUSIVE
-r                                                                            SUB_RERUNNABLE
-N                                                                            SUB_NOTIFY_END
-B                                                                            SUB_NOTIFY_BEGIN
-I                                                                            SUB_INTERACTIVE
-Ip                                                                           SUB_PTY
-Is                                                                           SUB_PTY_SHELL
-K                                                                            SUB2_BSUB_BLOCK*
-X "except_cond::action"          exceptList                                  SUB_EXCEPT
-T time_event                     timeEvent                                   SUB_TIME_EVENT
