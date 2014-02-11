package Genome::Logger;

use strict;
use warnings;

use Carp qw(croak);
use Log::Dispatch qw();
use Log::Dispatch::Screen::Color qw();
use Log::Dispatch::Screen qw();
use Memoize qw(memoize);

memoize('logger');
sub logger {
    assert_class_method(shift);

    my $logger = Log::Dispatch->new(@_);

    if (-t STDERR) {
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

sub screen {
    my $screen = Log::Dispatch::Screen->new(
        name => 'screen',
        min_level => 'info',
    );
    return $screen;
}

1;
