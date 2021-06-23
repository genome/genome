package Genome::Role::Logger;

use strict;
use warnings;

use Genome;
use UR::Role;

use Carp qw();
use Log::Dispatch qw();
use Log::Dispatch::Screen qw();
use Log::Dispatch::FileRotate qw();
use Module::Runtime qw(module_notional_filename use_package_optimistically);

my @log_levels;
BEGIN {
    @log_levels = keys %Log::Dispatch::CanonicalLevelNames;
    unless (@log_levels) {
        #older versions stored this differently
        @log_levels = keys %Log::Dispatch::LEVELS;
    }

    no warnings qw(redefine);
    sub Log::Dispatch::Base::_get_callbacks { return; }
}

role Genome::Role::Logger {
    has => [
        delegate_logger => {
            is => 'Log::Dispatch',
            is_constant => 1,
            is_calculated => 1,
            calculate => q( $self->_log_dispatch_init ),
        },
        screen => {
            is => 'Boolean',
            default => 1,
            doc => 'Display output to screen',
        },
        screen_level => {
            is => 'Text',
            default => 'warning',
            valid_values => \@log_levels,
            doc => 'The minimum level to display on the screen',
        },
        log_file => {
            is => 'FilePath',
            is_optional => 1,
            doc => 'Path to log file',
        },
        log_file_level => {
            is => 'Text',
            default => 'info',
            valid_values => \@log_levels,
            doc => 'Minimum level to display in the file log',
        },
        tie_stderr => {
            is => 'Boolean',
            default => 0,
            doc => '(warning) globally tie STDERR to this logger',
        },
    ],
};

sub _log_dispatch_init {
    my $self = shift;

    my $log = Log::Dispatch->new() || die "Can't create Log::Dispatch logger";

    if ($self->screen) {
        $log->add(_screen_to_add(min_level => $self->screen_level));
    }

    if ($self->log_file) {
        $log->add(
            Log::Dispatch::FileRotate->new(
                name => 'File',
                min_level => $self->log_file_level,
                filename => $self->log_file,
                mode => 'append',
                max => 5
            )
        );
    }

    if ($self->tie_stderr && !$self->screen) {
        tie(*STDERR, $self);
    } elsif ($self->tie_stderr && $self->screen) {
        $log->warning('WARNING: Disabling --tie-stderr since --screen is true.');
        $self->tie_stderr(0);
    }

    return $log;
}

# Create the logging functions info, notice, error, fatal, critical, etc.
BEGIN {
    for my $level ( @log_levels ) {
        my $name = join('::', __PACKAGE__, $level);
        my $log_sub = sub {
            my $self = shift;
            my $message = join(' ', @_);
            chomp $message;
            $message .= "\n";
            $self->delegate_logger->$level($message);
            return $message;
        };

        my $namef = $name . 'f';
        my $logf_sub = sub {
            my $self = shift;
            # sprintf inspects argument number
            my $message = sprintf(shift, @_);
            $self->$name($message);
        };

        no strict 'refs';
        *{$name} = $log_sub;
        *{$namef} = $logf_sub;
    }
}

sub stderror {
    my($self, $message) = @_;
    chomp $message;
    return $self->error("STDERR: $message\n");
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

# Thanks "socket puppet"
# http://stackoverflow.com/questions/8393667/is-there-a-way-to-redirect-prints-to-stderr-so-they-go-to-logdispatch

sub TIEHANDLE {
    my $logger = shift;
    die "$logger must use role Genome::Role::Logger" unless $logger->does('Genome::Role::Logger');
    return $logger;
}

sub PRINT {
    my $self = shift;
    $self->stderror(@_);
}

sub clear_logger {
    my $self = shift;
    $self->__invalidate_delegate_logger__;
}

sub _should_color_screen {
    return -t STDERR && has_color_screen_package();
}

sub has_color_screen_package {
    # Log::Dispatch::Screen::Color is not a "core" Log::Dispatch module
    my $name = use_package_optimistically('Log::Dispatch::Screen::Color');
    my $file = module_notional_filename($name);
    return $INC{$file};
}

sub _screen_to_add {
    if (_should_color_screen()) {
        return _create_color_screen(@_);
    } else {
        return _create_screen(@_);
    }
}

sub _create_screen {
    Log::Dispatch::Screen->new(
        name => 'screen',
        min_level => 'info',
        @_,
    );
}

sub _create_color_screen {
    Log::Dispatch::Screen::Color->new(
        name => 'screen',
        min_level => 'info',
        color => {
            debug => {
                text => 'magenta',
            },
            info => {}, # override default
            notice => {
                text => 'blue',
            },
            warning => {
                text => 'yellow',
            },
            error => {
                text => 'red',
            },
            critical => {
                text => 'red',
            },
            alert => {
                text => 'red',
            },
            emergency => {
                text => 'red',
            },
        },
        @_,
    );
}

BEGIN {
    for my $delegate_method (qw( add output remove ) ) {
        my $name = join('::', __PACKAGE__, $delegate_method);
        my $sub = sub {
            my $self = shift;
            $self->delegate_logger->$delegate_method(@_);
        };
        no strict 'refs';
        *{$name} = $sub;
    }
}

1;
