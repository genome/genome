package Genome::Utility::Instrumentation;
use parent 'Exporter';

# --- WARNING ---
# This is a prototype interface to the statsd server.
# It is likely to change, so use it at your own risk.
# --- WARNING ---


use strict;
use warnings;

our @EXPORT_OK = qw(
    dec
    decrement
    inc
    increment
    gauge
    timer
    timing
);

use Net::Statsd;
use Time::HiRes;

BEGIN {
    $Net::Statsd::HOST = $ENV{GENOME_STATSD_HOST} || '';
    $Net::Statsd::PORT = $ENV{GENOME_STATSD_PORT} || 0;
};


sub dec {
    eval {
        Net::Statsd::dec(@_);
    };
}

sub decrement {
    eval {
        Net::Statsd::decrement(@_);
    };
}


sub gauge {
    eval {
        Net::Statsd::gauge(@_);
    };
}


sub inc {
    eval {
        Net::Statsd::inc(@_);
    };
}

sub increment {
    eval {
        Net::Statsd::increment(@_);
    };
}


sub timer {
    my $code = pop @_;
    my @names = @_;

    my $start_time = Time::HiRes::time();

    eval {
        $code->();
    };
    if ($@) {
        my $error = $@;
        my $stop_time = Time::HiRes::time();
        for (@names) {
            timing("$_\_error", 1000 * ($stop_time - $start_time));
        }
        die $error;
    }

    my $stop_time = Time::HiRes::time();
    for(@names) {
        timing($_, 1000 * ($stop_time - $start_time));
    }
}

sub timing {
    eval {
        Net::Statsd::timing(@_);
    };
}

1;
