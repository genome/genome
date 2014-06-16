package Genome::Sys::Log;

use strict;
use warnings;

use Log::Log4perl qw(get_logger :levels);
use JSON;
use Digest::SHA qw(sha1);

# syslog, logstash, and UDP MTU can all affect this length.  Most recently we
# hit the UDP MTU limit going through the switches between blades and logstash.
# This is obviously site specific and should be move to configuration.
use constant MAX_SYSLOG_MSG_LEN => 1280;

my %MESSAGE_TYPE_TO_LOG_LEVEL = (
    debug_message => 'debug',
    status_message => 'info',
    warning_message => 'warn',
    error_message => 'error',
);

my $level = 0;
our %LOG_PRIORITY = map { $_ => ++$level } qw/debug info warn error fatal/;

# set with each callback fire for test cases
our $test_syslog;
our $last_syslog_level;
our $last_syslog_message;
our $last_syslog_retval;

my $log4perl;
my $callback = sub {
    my($self, $type, $message) = @_;

    # this should never happen given recent UR updates
    if (not defined $self) {
        no warnings;
        Carp::confess("self is undef, are you using the latest UR?: @_");
    }

    # make the logger on the first call
    unless ($log4perl) {
        $log4perl = get_logger();
        
        # get everything, since our real filter is on whether we choose to 
        # use the logger below.
        $log4perl->level($DEBUG);

        my $appender = Log::Log4perl::Appender->new(
            "Log::Dispatch::Syslog",
            ident => "GMS",
            facility => 'syslog',
        );
        $appender->layout(Log::Log4perl::Layout::PatternLayout->new('%m'));
        $log4perl->add_appender($appender);
    }

    my $level = $MESSAGE_TYPE_TO_LOG_LEVEL{$type};
    my $retval;

    # original logic runs only if the GENOME_SYS_LOG_DETAIL variable is not set

    unless ($ENV{GENOME_SYS_LOG_DETAIL}) {
        # by default we just log errors, and do so as text
        if ($level eq 'error') {
            if (substr($message,0,1) ne '{') {
                my $a = ref($self) ? $self->class . ' id('. $self->__display_name__.')' : $self;
                $message = $a . ': ' . $message;
            }

            # truncate message due to syslog limitation
            $message = substr($message, 0, MAX_SYSLOG_MSG_LEN);

            unless ($ENV{UR_DBI_NO_COMMIT}) {
                $retval = $log4perl->$level($message);
            }

            # these are used by test cases
            $last_syslog_level = $level;
            $last_syslog_message = $message;
            $last_syslog_retval = $retval;
        }
        return 1;
    }

    # detailed JSON logging occurs only when the GENOME_SYS_LOG_DETAIL variable is set for now
   
    # should we just standardize on JSON log entries?

    my $min_level = $ENV{GENOME_SYS_LOG_LEVEL};
    if (not $min_level) {
        # skip message b/c syslog level is not set at all
        return 1;
    }
    else {
        my $min_level_priority = Log::Log4perl::Level::to_priority(uc($min_level));
        my $this_level_priority = Log::Log4perl::Level::to_priority(uc($level));
        if ($this_level_priority < $min_level_priority) {
            return 1;
        }
    }

    my @stack_data;
    my @stack_txt;
    my $stack_level = 3;
    while(1) {
        my @c = caller($stack_level);
        if (@c == 0) {
            last;
        }
        else {
            my $module = $c[3];
            my $calling_line = $c[2];
            unshift @stack_data, \@c;
            unshift @stack_txt, "$calling_line|$module";
            $stack_level++;
        }
    }

    $stack_txt[-1] =~ s/\|.*//; 
    my $stack_txt = join(",", @stack_txt);
    my $stack_sha1 = sha1($stack_txt);

    my $offset = ($stack_data[-1][0] eq 'Genome::Sys' ? -2 : -1);
    my $package = $stack_data[$offset][0];
    my $filename = $stack_data[$offset][1];
    my $line = $stack_data[$offset][2];
    my $subroutine = ($stack_data[$offset-1] ? $stack_data[$offset-1][3] : '');

    # if the message is json simply extend it
    my $incoming_json_data;
    if (substr($message,0,1) eq '{') {
        $incoming_json_data = decode_json $message;
    }

    # use brief keys since we are space-constrained in syslog
    my $msg = do {
        no warnings;
        {
            f => $filename,
            l => $line,
            m => $subroutine,
            c => $stack_txt,
            s => $stack_sha1,
            h => Sys::Hostname::hostname(),
            p => $$,
            j => $ENV{LSB_JOBID},
            u => $ENV{USER},
            b => $ENV{GENOME_BUILD_ID},
            type => $level,
            ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} ? (test => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME}) : ()),
            ($incoming_json_data ? %$incoming_json_data : (msg => $message)),
        };
    };
    my $json = encode_json $msg;

    if (length($json) > MAX_SYSLOG_MSG_LEN) {
        # omit the explicit stack if the message is too long
        # the sha1 of the stack is still there
        delete $msg->{c};
        $json = encode_json $msg;
    }
    
    if ($test_syslog or $ENV{UR_DBI_NO_COMMIT}) { 
        # remove this stderr printing if we make JSON logging default
        print STDERR "LOGGING ($level) (test/skip): $json\n";
        # TODO: it would be better to redirect the appender instead of skipping using it.
        $retval = 1;
    }
    else {
        # remove this stderr printing if we make JSON logging default
        print STDERR "LOGGING ($level): $json\n";
        $retval = $log4perl->$level($json);
    }
    
    # these are used by test cases
    $last_syslog_level = $level;
    $last_syslog_message = $json;
    $last_syslog_retval = $retval;

    return $retval;
};

for my $type (sort keys %MESSAGE_TYPE_TO_LOG_LEVEL) {
    UR::Object->add_observer(
        aspect => $type,
        callback => $callback, 
    );
}

1;

