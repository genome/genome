package Genome::Sys::LSF::bsub;

use strict;
use warnings;

use Exporter qw(import);
use IPC::System::Simple qw(capture);
use Params::Validate qw(:types);

our @EXPORT = qw(bsub);
our @EXPORT_OK = qw(bsub);

sub bsub {
    my @command = _command_builder(@_);

    # lazy load so we don't break /gsc/bin/perl (until we have to)
    require IPC::System::Simple;
    my @output = IPC::System::Simple::capture(@command);

    my $job_id = ($output[-1] =~ /^Job <(\d+)> is submitted to/)[0];
    unless ($job_id) {
        die "Could not get job id from bsub output!";
    }

    return $job_id;
}

sub _command_builder {
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
        push @args, map { _option_mapper($_) => $args{$_} } sort keys %args;
    }

    if (ref($command) ne 'ARRAY') {
        $command = [$command];
    }

    return ('bsub', @args, @$command);
}

sub _option_mapper {
    my $arg = shift;
    return {
        email           => '-u',
        err_file        => '-e',
        hold_job        => '-H',
        job_group       => '-g',
        log_file        => '-o',
        project         => '-P',
        queue           => '-q',
        send_job_report => '-N',
    }->{$arg};
}

sub _args {
    my %args = Params::Validate::validate(
        @_, {
            _flag_spec(),
            queue => {
                optional => 1,
                type => SCALAR,
                callbacks => {
                    'valid LSF queue' => \&_valid_lsf_queue,
                },
            },
            project => {
                optional => 1,
                type => SCALAR,
            },
            email => {
                optional => 1,
                type => SCALAR,
            },
            err_file => {
                optional => 1,
                type => SCALAR,
            },
            log_file => {
                optional => 1,
                type => SCALAR,
            },
            job_group => {
                optional => 1,
                type => SCALAR,
            },
            cmd => {
                type => SCALAR | ARRAYREF,
            },
        }
    );
    return %args;
}

sub _flags {
    return sort qw(
        hold_job
        send_job_report
    );
}

sub _flag_spec {
    return map {
        $_ => {
            optional => 1,
            type => BOOLEAN,
        }
    } _flags();
}

sub _valid_lsf_queue {
    return grep { /$_[0]/ } _queues();
}

sub _queues {
    my @output = capture('bqueues', '-l');
    my @queues = map { (/^QUEUE:\s+(\S+)/)[0] } @output;
    return @queues;
}

1;
