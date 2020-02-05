package Genome::Sys::LSF::bsub;

use strict;
use warnings;

use Genome::Sys;
use Exporter qw(import);
use Params::Validate qw(:types);
use List::MoreUtils qw(any);
use Genome::Sys::LSF::ResourceParser qw();

our @EXPORT = qw(bsub);
our @EXPORT_OK = qw(bsub);

sub run {
    my $executable = shift;
    my @args = args_builder(@_);

    if (ref($executable) ne 'ARRAY') {
        $executable = [$executable];
    }

    local $ENV{LSB_SUB_ADDITIONAL} = Genome::Config::get('lsb_sub_additional') || $ENV{LSB_SUB_ADDITIONAL};
    local $ENV{LSF_DOCKER_VOLUMES} = Genome::Config::get('docker_volumes') || $ENV{LSF_DOCKER_VOLUMES};

    my @output = Genome::Sys->capture(@$executable, @args);

    my $job_id = ($output[-1] =~ /^Job <(\d+)> is submitted to/)[0];
    unless ($job_id) {
        die "Could not get job id from bsub output!";
    }

    return $job_id;
}

sub bsub {
    return run('bsub', @_);
}

sub args_builder {
    my %args = _args(@_);
    my $command = delete $args{cmd};
    my $resource_string = delete $args{resource_string};

    my @args;
    for my $flag (_flags()) {
        if ($args{$flag}) {
            push @args, _option_mapper($flag);
        }
        delete $args{$flag};
    }
    if (keys %args) {
        push @args, map { _option_mapper($_), $args{$_} } sort keys %args;
    }

    if (ref($command) ne 'ARRAY') {
        $command = [$command];
    }

    if ($resource_string) {
        my $parsed = Genome::Sys::LSF::ResourceParser::parse_lsf_params($resource_string);
        if (exists $parsed->{options}{resReq}) {
            push @args, '-R', $parsed->{options}{resReq};
        }
        if (exists $parsed->{options}{numProcessors}) {
            push @args, '-n', $parsed->{options}{numProcessors};
        }
        if (exists $parsed->{rLimits}{RSS}) {
            push @args, '-M', $parsed->{rLimits}{RSS};
        }
    }

    return (@args, @$command);
}

sub _args {
    return Params::Validate::validate(@_, _args_spec());
}

sub _option_mapper {
    my $arg = shift;
    my $spec = _args_spec();
    if (not exists $spec->{$arg}{option_flag}) {
        die qq(Could not find option flag for '$arg' in args spec);
    }
    return $spec->{$arg}{option_flag};
}

sub _flags {
    my $spec = _args_spec();
    return grep { $spec->{$_}{type} == BOOLEAN } sort keys %{$spec};
}

sub _args_spec {
    return {
        rerunnable => {
            optional => 1,
            type => BOOLEAN,
            option_flag => '-r',
        },
        never_rerunnable => {
            optional => 1,
            type => BOOLEAN,
            option_flag => '-rn',
        },
        hold_job => {
            optional => 1,
            option_flag => '-H',
            type => BOOLEAN,
        },
        send_job_report => {
            optional => 1,
            option_flag => '-N',
            type => BOOLEAN,
        },
        queue => {
            optional => 1,
            option_flag => '-q',
            type => SCALAR,
            callbacks => {
                'valid LSF queue' => \&_valid_lsf_queue,
            },
        },
        project => {
            optional => 1,
            option_flag => '-P',
            type => SCALAR,
        },
        email => {
            optional => 1,
            option_flag => '-u',
            type => SCALAR,
        },
        err_file => {
            optional => 1,
            option_flag => '-e',
            type => SCALAR,
        },
        log_file => {
            optional => 1,
            option_flag => '-o',
            type => SCALAR,
        },
        job_group => {
            optional => 1,
            option_flag => '-g',
            type => SCALAR,
        },
        job_name => {
            optional => 1,
            option_flag => '-J',
            type => SCALAR,
        },
        depend_on => {
            optional => 1,
            option_flag => '-w',
            type => SCALAR,
        },
        wait_for_completion => {
            optional => 1,
            option_flag => '-K',
            type => BOOLEAN,
        },
        post_exec_cmd => {
            optional => 1,
            option_flag => '-Ep',
            type => SCALAR,
        },
        cmd => {
            type => SCALAR | ARRAYREF,
        },
        resource_string => {
            optional => 1,
            type => SCALAR,
        },
    };
}

sub _valid_lsf_queue {
    my $requested_queue = shift;

    system(qq(bqueues $requested_queue 1> /dev/null 2>&1));
    return $? == 0;
}

1;
