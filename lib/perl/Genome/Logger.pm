package Genome::Logger;

use strict;
use warnings;

use Carp qw(croak);
use Log::Dispatch qw();
use Log::Dispatch::Screen qw();
use Memoize qw(memoize);
use Module::Runtime qw(module_notional_filename use_package_optimistically);

memoize('logger');
sub logger {
    assert_class_method(shift);

    my $logger = Log::Dispatch->new(@_);

    if (should_color_screen()) {
        $logger->add(color_screen());
    } else {
        $logger->add(screen());
    }

    return $logger;
}

sub assert_class_method_error { 'Must be called as class method' }
sub assert_class_method {
    my $class = shift;

    # to ensure memoize works we are strict about this
    unless ($class && $class eq __PACKAGE__) {
        croak assert_class_method_error();
    }
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
    );
}

my @levels = keys %Log::Dispatch::LEVELS;
for my $level (@levels) {
    my $name = join('::', __PACKAGE__, $level);
    no strict 'refs';
    *{$name} = sub {
        my $class = shift;
        return $class->logger->$level(@_);
    };
}

sub screen {
    my $screen = Log::Dispatch::Screen->new(
        name => 'screen',
        min_level => 'info',
    );
    return $screen;
}

1;
