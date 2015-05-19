package Genome::Role::Logger;

use strict;
use warnings;

use Genome;
use Log::Dispatch;
use Log::Dispatch::Screen;
use Log::Dispatch::FileRotate;

my @log_levels = keys %Log::Dispatch::LEVELS;

class Genome::Role::Logger {
    has => [
        screen => {
            is => 'Boolean',
            default => 1,
            doc => 'Display output to screen.',
        },
        screen_level => {
            is => 'Text',
            default => 'warning',
            valid_values => \@log_levels,
            doc => 'The minimum level to display on the scren.',
        },
        log_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'Path to log file.',
        },
        log_file_level => {
            is => 'Text',
            default => 'info',
            valid_values => \@log_levels,
            doc => 'The minimum level to display in the log.',
        },
        tie_stderr => {
            is => 'Boolean',
            default => 0,
            doc => '(warning) globally tie STDERR to this logger',
        },
        log_dispatch => {
            is => 'Log::Dispatch',
            is_calculated => 1,
            is_constant => 1,
            calculate => q($self->log_dispatch_init),
        },
    ],
};

sub log_dispatch_init {
    my $self = shift;

    my $log = Log::Dispatch->new() || die;

    if ($self->screen) {
        $log->add(Genome::Logger->screen_to_add(min_level => $self->screen_level));
    }

    if ($self->log_file) {
        $log->add(
            Log::Dispatch::FileRotate->new(
                name => 'File',
                min_level => $self->log_file_level,
                filename => $self->log_file,
                mode => 'append',
                max => 5,
            )
        );
    }

    if ($self->tie_stderr && !$self->screen) {
        tie(*STDERR, $self);
    } elsif ($self->tie_stderr && $self->screen) {
        $log->warning("WARNING: Disabling --tie_stderr since --screen is true.\n");
        $self->tie_stderr(0);
    }

    return $log;
}

##############################################
# This is almost duplicated in Genome::Logger.
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# create object methods for each log level
for my $log_level (@log_levels) {
    my $name = join('::', __PACKAGE__, $log_level);
    my $namef = $name . 'f';
    no strict 'refs';
    *{$name} = sub {
        my $self = shift;
        return $self->log_dispatch->$log_level(@_);
    };
    *{$namef} = sub {
        my $self = shift;
        # sprintf inspects argument number
        my $message = sprintf(shift, @_);
        $self->$name($message);
    };
}

sub croak {
    my $self = shift;
    my $level = shift;

    unless ($self->can($level)) {
        Carp::croak "invalid level: $level";
    }

    Carp::croak $self->$level(@_);
}

sub fatal {
    my $self = shift;
    $self->croak('critical', @_);
}

sub fatalf {
    my $self = shift;
    $self->croak('criticalf', @_);
}

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# End almost duplication from Genome::Logger.
#############################################

sub stderror {
    my ($self, $message) = @_;
    chomp $message;
    $message = uc('stderr') . ": $message\n";
    return $self->log_dispatch->error($message);
}

# Thanks "socket puppet"
# http://stackoverflow.com/questions/8393667/is-there-a-way-to-redirect-prints-to-stderr-so-they-go-to-logdispatch

sub TIEHANDLE {
    my $logger = shift;
    die unless $logger->isa('Genome::Role::Logger');
    return $logger;
}

sub PRINT {
    my $self = shift;
    $self->stderror(@_);
}

1;
