package Genome::Sys::LSF::ResourceParser;

use strict;
use warnings FATAL => qw(all);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(parse_lsf_params parse_resource_requirements);

use Getopt::Long qw(GetOptionsFromString);
use IO::String;

sub parse_resource_requirements {
    my ($res_req) = @_;
    $res_req = '' unless defined $res_req;

    # Platform LSF Admin Guide
    # About resource requirement strings

    # Simple syntax
    # select[selection_string] order[order_string] rusage[usage_string [, usage_string]
    # [|| usage_string] ...] span[span_string] same[same_string] cu[cu_string] affinity[affinity_string]
    my $simple_string = qr/
        \s*
            \w+\[ [^\]]+ \]
        (\s+\w+\[ [^\]]+ \])*
        \s*
    /xms;

    # Compound syntax
    # num1*{simple_string1} + num2*{simple_string2} + ...
    my $compound_string = qr/
        \d+\*{$simple_string}
        (\s*\+\s*\d+\*{$simple_string})*
    /xms;

    my $admissable_res_req = qr/^( | $simple_string | $compound_string )$/xms;
    return $res_req =~ $admissable_res_req;
}

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

    my $res_req = exists($lsf_params{options}{resReq})
        ? $lsf_params{options}{resReq} : '';

    if ($parse_ok && !parse_resource_requirements($res_req)) {
        return (0, \%lsf_params, 'Invalid resource requirements specification');
    }
    else {
        return ($parse_ok, \%lsf_params, $message);
    }
}

sub _get_options_from_string {
    my $getopt_fh = IO::String->new;
    local *STDERR = $getopt_fh;
    my ($ret, $args) = GetOptionsFromString(@_);
    my $parse_ok = $ret && (scalar @$args == 0);
    return ($parse_ok, ${$getopt_fh->string_ref});
}

sub _create_getopt_specs {
    my ($destination, %option) = @_;

    my @getopt_specs;
    while (my ($bsub_option, $submit_field) = each %option) {
        push @getopt_specs, _lookup_spec_for_bsub_option(
            $bsub_option, $submit_field, $destination);
    }

    return @getopt_specs;
}

sub _lookup_spec_for_bsub_option {
    my ($bsub_option, $submit_field, $destination) = @_;
    my $formatter = _formatters($bsub_option);
    if ($bsub_option eq 'n') {
        return("$bsub_option=s" => sub {
                my ($option_name, $option_value) = @_;
                my ($num, $max_num) = split q{,}, $option_value;
                $destination->{numProcessors} = $formatter->($num);
                $destination->{maxNumProcessors} = $formatter->($max_num)
                    if defined $max_num;
            });
    }
    else {
        return("$bsub_option=s" => sub {
                my ($option_name, $option_value) = @_;
                $destination->{$submit_field} = $formatter->($option_value);
            });
    }
}

sub _create_option_spec {
    my ($option, $dictionary) = @_;
    return ("$option=s" => sub { my ($option_name, $value, )});
}

sub string_formatter {
    my $value = shift;
    return "$value";
}

sub number_formatter {
    my $value = shift;
    return $value + 0;
}

sub minutes_to_seconds_formatter {
    my $value = shift;
    return $value * 60;
}

my %formatters = (
    'b' => \&string_formatter,
    'e' => \&string_formatter,
    'g' => \&string_formatter,
    'i' => \&string_formatter,
    'J' => \&string_formatter,
    'u' => \&string_formatter,
    'n' => \&number_formatter,
    'o' => \&string_formatter,
    'E' => \&string_formatter,
    'Ep' => \&string_formatter,
    'P' => \&string_formatter,
    'q' => \&string_formatter,
    'R' => \&string_formatter,
    't' => \&string_formatter,

    'c' => \&minutes_to_seconds_formatter,
    'M' => \&number_formatter,
    'F' => \&number_formatter,
    'S' => \&number_formatter,
);
sub _formatters {
    my $key = shift;
    return $formatters{$key};
}

my %valid_options = (
    'b' => 'beginTime',
    'e' => 'errFile',
    'g' => 'jobGroup',
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
