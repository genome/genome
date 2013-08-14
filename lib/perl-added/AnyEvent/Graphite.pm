package AnyEvent::Graphite;

use warnings;
use strict;

our $VERSION = '0.08';

use AnyEvent;
use AnyEvent::Socket;
use AnyEvent::Handle;

sub new {
    my($class, %args) = @_;

    # the rest of the args will get passed to the object
    # so remove host and port as we'll possibly set them manually here
    my $host = delete $args{host} || '127.0.0.1';
    my $port = delete $args{port} || 2003;

    bless {
        host => $host,
        port => $port,
        %args,
    }, $class;
}

sub send {
    my($self, $id, $value, $ts) = @_;
    return unless (defined($id) && defined($value));
    $ts ||= AE::now; # cached time() from the framework
    if($self->{conn}) {
        $self->{conn}->push_write(join(" ", $id, $value, $ts) . "\n");
    } else {
        my $handle; $handle = new AnyEvent::Handle
            connect => [$self->{host} => $self->{port}],
            on_error => sub {
                # we have no reason to exist w/out a connection.
                die "Unable to connect to Graphite server at $self->{host}:$self->{port}: $!\n";
            };
        $self->{conn} = $handle;
        $self->send($id, $value, $ts);
    }
}

sub finish {
    # you need to let the client actually process the data, so the event loop needs to run
    # this lets you let your data get out, but still exit when it's done
    my($self) = @_;
    $self->{conn}->on_drain(sub {
        print "All data transmitted\n";
        exit;
    });
    AnyEvent->condvar->recv;
}

=head1 NAME

AnyEvent::Graphite - A non-blocking Graphite (http://graphite.wikidot.com/) client

=head1 VERSION

Version 0.08

=head1 SYNOPSIS

    my $graphite = AnyEvent::Graphite->new(
        host => '127.0.0.1',
        port => '2003',
    );
    $graphite->send("a.b.c.d", $value, $timestamp); #timestamp is optional

=head1 AUTHOR

Joshua Barratt, C<< <josh at mediatemple.net> >>

=head1 COPYRIGHT & LICENSE

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of AnyEvent::Graphite
