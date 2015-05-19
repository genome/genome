package Genome::Logger;

use strict;
use warnings;

use Carp qw();
use Log::Dispatch qw();
use Log::Dispatch::Screen qw();
use Module::Runtime qw(module_notional_filename use_package_optimistically);

my $logger;
sub logger {
    my $class = shift;
    if ($logger) {
        return $logger;
    }

    $logger = Log::Dispatch->new(@_);
    $logger->add(screen_to_add());

    return $logger;
}
sub clear_logger {
    $logger = undef;
}

sub should_color_screen {
    return -t STDERR && has_color_screen_package();
}

sub has_color_screen_package {
    # Log::Dispatch::Screen::Color is not a "core" Log::Dispatch module
    my $name = use_package_optimistically('Log::Dispatch::Screen::Color');
    my $file = module_notional_filename($name);
    return $INC{$file};
}

sub screen_to_add {
    my $class = shift;
    if (should_color_screen()) {
        return color_screen(@_);
    } else {
        return screen(@_);
    }
}

sub screen {
    my $screen = Log::Dispatch::Screen->new(
        name => 'screen',
        min_level => 'info',
        @_,
    );
    return $screen;
}

sub color_screen {
    my $screen = Log::Dispatch::Screen::Color->new(
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

####################################################
# This is almost duplicated in Genome::Role::Logger.
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

my @levels = keys %Log::Dispatch::LEVELS;
for my $level (@levels) {
    my $name = join('::', __PACKAGE__, $level);
    my $namef = $name . 'f';
    no strict 'refs';
    *{$name} = sub {
        my $class = shift;
        $class->logger->$level(@_);
        return join(' ', @_);
    };
    *{$namef} = sub {
        my $class = shift;
        # sprintf inspects argument number
        my $message = sprintf(shift, @_);
        $class->$name($message);
    };
}

sub croak {
    my $class = shift;
    my $level = shift;

    unless ($class->can($level)) {
        Carp::croak "invalid level: $level";
    }

    Carp::croak $class->$level(@_);
}

sub fatal {
    my $class = shift;
    $class->croak('critical', @_);
}

sub fatalf {
    my $class = shift;
    $class->croak('criticalf', @_);
}

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# End almost duplication from Genome::Role::Logger.
###################################################

1;
