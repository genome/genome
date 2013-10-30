package Genome::Utility::Instrumentation;


# --- WARNING ---
# This is a prototype interface to the statsd server.
# It is likely to change, so use it at your own risk.
# --- WARNING ---


use strict;
use warnings;

use Genome;

use Net::Statsd;
use Time::HiRes;

BEGIN {
    $Net::Statsd::HOST = $ENV{GENOME_STATSD_HOST} || '';
    $Net::Statsd::PORT = $ENV{GENOME_STATSD_PORT} || 0;
};


sub dec {
    return unless $Net::Statsd::HOST;
    eval {
        Net::Statsd::dec(@_);
    };
}

sub decrement {
    return unless $Net::Statsd::HOST;
    eval {
        Net::Statsd::decrement(@_);
    };
}


sub gauge {
    return unless $Net::Statsd::HOST;
    eval {
        Net::Statsd::gauge(@_);
    };
}


sub inc {
    return unless $Net::Statsd::HOST;
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
    return unless $Net::Statsd::HOST;
    my ($name, $code) = @_;

    my $start_time = Time::HiRes::time();

    eval {
        $code->();
    };
    if ($@) {
        my $error = $@;
        my $stop_time = Time::HiRes::time();
        eval {
            my $error_name = "$name\_error";
            Net::Statsd::timing($error_name, 1000 * ($stop_time - $start_time));
        };
        die $error;
    }

    my $stop_time = Time::HiRes::time();
    eval {
        Net::Statsd::timing($name, 1000 * ($stop_time - $start_time));
    };
}

sub timing {
    return unless $Net::Statsd::HOST;
    eval {
        Net::Statsd::timing(@_);
    };
}

1;
