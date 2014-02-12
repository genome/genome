package Genome::Sys::LSF::bsub;

use strict;
use warnings;

use Exporter qw(import);
use Params::Validate qw(:types);

our @EXPORT = qw(bsub);
our @EXPORT_OK = qw(bsub);

sub run {
    my $executable = shift;
    my @args = args_builder(@_);

    if (ref($executable) ne 'ARRAY') {
        $executable = [$executable];
    }

    my @output = _capture(@$executable, @args);

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
        push @args, map { _option_mapper($_) => $args{$_} } sort keys %args;
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
        cmd => {
            type => SCALAR | ARRAYREF,
        },
    };
}

sub _valid_lsf_queue {
    return grep { /$_[0]/ } _queues();
}

sub _queues {
    my @output = _capture('bqueues', '-l');
    my @queues = map { (/^QUEUE:\s+(\S+)/)[0] } @output;
    return @queues;
}

sub _capture {
    # lazy load so we don't break /gsc/bin/perl (until we have to)
    require IPC::System::Simple;
    return IPC::System::Simple::capture(@_);
}

1;
