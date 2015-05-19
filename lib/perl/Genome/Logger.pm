package Genome::Logger;

use strict;
use warnings;

use Carp qw();
use Log::Dispatch qw();
use Log::Dispatch::Screen qw();
use Module::Runtime qw(module_notional_filename use_package_optimistically);

use UR;

class Genome::Logger {
    has => {
        delegate_logger => {
            is => 'Log::Dispatch',
            is_constant => 1,
        },
    },
};

my $logger;
sub logger {
    my $class = shift;
    if ($logger) {
        return $logger;
    }

    $logger = Genome::Logger->create(
        delegate_logger => Log::Dispatch->new(@_),
    );
    $logger->delegate_logger->add(screen_to_add());

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

sub normalize_self {
    my $class = shift;
    my $self = ref $class ? $class : $class->logger;
    return $self;
}

my @levels = keys %Log::Dispatch::LEVELS;
for my $level (@levels) {
    my $name = join('::', __PACKAGE__, $level);
    my $namef = $name . 'f';
    no strict 'refs';
    *{$name} = sub {
        my $self = normalize_self(shift);
        $self->delegate_logger->$level(@_);
        return join(' ', @_);
    };
    *{$namef} = sub {
        my $self = normalize_self(shift);
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

for my $delegate_method (qw(add output remove)) {
    my $name = join('::', __PACKAGE__, $delegate_method);
    no strict 'refs';
    *{$name} = sub {
        my $self = shift;
        return $self->delegate_logger->$delegate_method(@_);
    }
}

1;
