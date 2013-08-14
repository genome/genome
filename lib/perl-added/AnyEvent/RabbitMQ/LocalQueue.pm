package AnyEvent::RabbitMQ::LocalQueue;

use strict;
use warnings;

our $VERSION = '1.08';

sub new {
    my $class = shift;
    return bless {
        _message_queue    => [],
        _drain_code_queue => [],
    }, $class;
}

sub push {
    my $self = shift;

    CORE::push @{$self->{_message_queue}}, @_;
    return $self->_drain_queue();
}

sub get {
    my $self = shift;

    CORE::push @{$self->{_drain_code_queue}}, @_;
    return $self->_drain_queue();
}

sub _drain_queue {
    my $self = shift;

    my $message_count = scalar @{$self->{_message_queue}};
    my $drain_code_count = scalar @{$self->{_drain_code_queue}};

    my $count = $message_count < $drain_code_count
              ? $message_count : $drain_code_count;

    for (1 .. $count) {
        &{shift @{$self->{_drain_code_queue}}}(
            shift @{$self->{_message_queue}}
        );
    }

    return $self;
}

sub _flush {
    my ($self, $frame) = @_;

    $self->_drain_queue;

    while (my $cb = shift @{$self->{_drain_code_queue}}) {
        local $@; # Flush frames immediately, throwing away errors for on-close
        eval { $cb->($frame) };
    }
}

1;

