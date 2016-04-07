package Genome::Sys::LSF::bsub;

use strict;
use warnings;

use Genome::Sys;
use Exporter qw(import);
use Params::Validate qw(:types);
use List::MoreUtils qw(any);
use Genome::Utility::Email;

our @EXPORT = qw(bsub);
our @EXPORT_OK = qw(bsub);

sub run {
    my $executable = shift;
    my @args = args_builder(@_);

    if (ref($executable) ne 'ARRAY') {
        $executable = [$executable];
    }

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
    };
}

sub _valid_lsf_queue {
    my $requested_queue = shift;

    my $username = Genome::Sys->username;
    if ($username eq 'apipe-builder' and $requested_queue eq 'apipe') {
        my $message = join("\n",
                            'apipe-builder using apipe queue',
                            'Host ' . $ENV{HOSTNAME},
                            'LSF jobID ' . $ENV{LSB_JOBID},
                            'submission host '. $ENV{LSB_SUB_HOST});

        Genome::Utility::Email::send(
            from => 'abrummet@genome.wustl.edu',
            to => 'abrummet@genome.wustl.edu',
            subject => 'apipe-builder using apipe queue',
            body => Carp::longmess($message),
        );
    }
    return any { $requested_queue eq $_ } _queues();
}

sub _queues {
    my @output = Genome::Sys->capture('bqueues', '-l');
    my @queues = map { (/^QUEUE:\s+(\S+)/)[0] } @output;
    return @queues;
}

1;
