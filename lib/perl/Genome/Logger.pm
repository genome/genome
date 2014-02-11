package Genome::Logger;

use strict;
use warnings;

use Carp qw(croak);
use Log::Dispatch qw();
use Log::Dispatch::Screen::Color qw();
use Memoize qw(memoize);

memoize('logger');
sub logger {
    assert_class_method(shift);

    my $logger = Log::Dispatch->new(@_);
    $logger->add(color_screen());
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

1;
