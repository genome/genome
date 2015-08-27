package Genome::Logger;

use strict;
use warnings;

use Genome;

class Genome::Logger {
    is => 'UR::Singleton',
    roles => 'Genome::Role::Logger',
};

my $logger;
sub delegate_logger : Overrides(Genome::Role::Logger) {
    my $self = shift;
    unless ($logger) {
        my $self = $self->_singleton_object;
        $logger = $self->_log_dispatch_init();
    }
    return $logger;
}

*logger = \&get;

sub clear_logger : Overrides(Genome::Role::Logger) {
    undef $logger;
}

1;
