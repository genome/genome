package AnyEvent::RabbitMQ::Channel;

use strict;
use warnings;

use Scalar::Util qw(weaken);
use AnyEvent::RabbitMQ::LocalQueue;

our $VERSION = '1.08';

sub new {
    my $class = shift;
    my $self = bless {
        @_, # id, connection, on_close
        _is_open       => 0,
        _is_active     => 0,
        _queue         => AnyEvent::RabbitMQ::LocalQueue->new,
        _content_queue => AnyEvent::RabbitMQ::LocalQueue->new,
        _consumer_cbs  => {},
        _return_cbs    => {},
    }, $class;
    weaken($self->{connection});
    return $self;
}

sub queue {
    my $self = shift;
    return $self->{_queue};
}

sub open {
    my $self = shift;
    my %args = @_;

    if ($self->{_is_open}) {
        $args{on_failure}->('Channel has already been opened');
        return $self;
    }

    $self->{connection}->_push_write_and_read(
        'Channel::Open', {}, 'Channel::OpenOk',
        sub {
            $self->{_is_open}   = 1;
            $self->{_is_active} = 1;
            $args{on_success}->();
        },
        sub {
            $args{on_failure}->(@_);
        },
        $self->{id},
    );

    return $self;
}

sub close {
    my $self = shift;
    my $connection = $self->{connection}
        or return;
    my %args = $connection->_set_cbs(@_);

    # Ensure to remove this channel from the connection even if we're not
    # fully open to ensure $rf->close works always.
    # FIXME - We can end up racing here so the server thinks the channel is
    # open, but we've closed it - a more elegant fix would be to mark that
    # the channel is opening, and wait for it to open before closing it
    if (!$self->{_is_open}) {
        $self->{connection}->delete_channel($self->{id});
        $args{on_success}->($self);
        return $self;
    }

    return $self->_close(%args) if 0 == scalar keys %{$self->{_consumer_cbs}};

    for my $consumer_tag (keys %{$self->{_consumer_cbs}}) {
        $self->cancel(
            consumer_tag => $consumer_tag,
            on_success   => sub {
                $self->_close(%args);
            },
            on_failure   => sub {
                $self->_close(%args);
                $args{on_failure}->(@_);
            }
        );
    }

    return $self;
}

sub _close {
    my $self = shift;
    my %args = @_;

    $self->{connection}->_push_write_and_read(
        'Channel::Close', {}, 'Channel::CloseOk',
        sub {
            $self->{_is_open} = 0;
            $self->{_is_active} = 0;
            $self->{connection}->delete_channel($self->{id});
            $args{on_success}->();
        },
        sub {
            $self->{_is_open} = 0;
            $self->{_is_active} = 0;
            $self->{connection}->delete_channel($self->{id});
            $args{on_failure}->();
        },
        $self->{id},
    );

    return $self;
}

sub declare_exchange {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Exchange::Declare',
        {
            type        => 'direct',
            passive     => 0,
            durable     => 0,
            auto_delete => 0,
            internal    => 0,
            %args, # exchange
            ticket      => 0,
            nowait      => 0, # FIXME
        },
        'Exchange::DeclareOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub delete_exchange {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Exchange::Delete',
        {
            if_unused => 0,
            %args, # exchange
            ticket    => 0,
            nowait    => 0, # FIXME
        },
        'Exchange::DeleteOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub declare_queue {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Queue::Declare',
        {
            queue       => '',
            passive     => 0,
            durable     => 0,
            exclusive   => 0,
            auto_delete => 0,
            no_ack      => 1,
            %args,
            ticket      => 0,
            nowait      => 0, # FIXME
        },
        'Queue::DeclareOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );
}

sub bind_queue {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Queue::Bind',
        {
            %args, # queue, exchange, routing_key
            ticket => 0,
            nowait => 0, # FIXME
        },
        'Queue::BindOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub unbind_queue {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Queue::Unbind',
        {
            %args, # queue, exchange, routing_key
            ticket => 0,
        },
        'Queue::UnbindOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub purge_queue {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Queue::Purge',
        {
            %args, # queue
            ticket => 0,
            nowait => 0, # FIXME
        },
        'Queue::PurgeOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub delete_queue {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Queue::Delete',
        {
            if_unused => 0,
            if_empty  => 0,
            %args, # queue
            ticket    => 0,
            nowait    => 0, # FIXME
        },
        'Queue::DeleteOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub publish {
    my $self = shift;
    my %args = @_;

    return $self if !$self->{_is_active};

    my $header_args = delete $args{header}    || {};
    my $body        = delete $args{body}      || '';
    my $return_cb   = delete $args{on_return} || sub {};

    $self->_publish(
        %args,
    )->_header(
        $header_args, $body,
    )->_body(
        $body,
    );

    return $self if !$args{mandatory} && !$args{immediate};

    $self->{_return_cbs}->{
        ($args{exchange} || '') . '_' . $args{routing_key}
    } = $return_cb;

    return $self;
}

sub _publish {
    my $self = shift;
    my %args = @_;

    $self->{connection}->_push_write(
        Net::AMQP::Protocol::Basic::Publish->new(
            exchange  => '',
            mandatory => 0,
            immediate => 0,
            %args, # routing_key
            ticket    => 0,
        ),
        $self->{id},
    );

    return $self;
}

sub _header {
    my ($self, $args, $body,) = @_;

    $self->{connection}->_push_write(
        Net::AMQP::Frame::Header->new(
            weight       => $args->{weight} || 0,
            body_size    => length($body),
            header_frame => Net::AMQP::Protocol::Basic::ContentHeader->new(
                content_type     => 'application/octet-stream',
                content_encoding => undef,
                headers          => {},
                delivery_mode    => 1,
                priority         => 1,
                correlation_id   => undef,
                expiration       => undef,
                message_id       => undef,
                timestamp        => time,
                type             => undef,
                user_id          => $self->{connection}->login_user,
                app_id           => undef,
                cluster_id       => undef,
                %$args,
            ),
        ),
        $self->{id},
    );

    return $self;
}

sub _body {
    my ($self, $body,) = @_;

    $self->{connection}->_push_write(
        Net::AMQP::Frame::Body->new(payload => $body),
        $self->{id},
    );

    return $self;
}

sub consume {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    my $consumer_cb = delete $args{on_consume} || sub {};
    
    $self->{connection}->_push_write_and_read(
        'Basic::Consume',
        {
            consumer_tag => '',
            no_local     => 0,
            no_ack       => 1,
            exclusive    => 0,
            %args, # queue
            ticket       => 0,
            nowait       => 0, # FIXME
        },
        'Basic::ConsumeOk', 
        sub {
            my $frame = shift;
            $self->{_consumer_cbs}->{
                $frame->method_frame->consumer_tag
            } = $consumer_cb;
            $cb->($frame);
        },
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub cancel {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    if (!defined $args{consumer_tag}) {
        $failure_cb->('consumer_tag is not set');
        return $self;
    }

    if (!$self->{_consumer_cbs}->{$args{consumer_tag}}) {
        $failure_cb->('Unknown consumer_tag');
        return $self;
    }

    $self->{connection}->_push_write_and_read(
        'Basic::Cancel',
        {
            %args, # consumer_tag
            nowait => 0,
        },
        'Basic::CancelOk', 
        sub {
            my $frame = shift;
            delete $self->{_consumer_cbs}->{$args{consumer_tag}};
            $cb->($frame);
        },
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub get {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Basic::Get',
        {
            no_ack => 1,
            %args, # queue
            ticket => 0,
        },
        [qw(Basic::GetOk Basic::GetEmpty)], 
        sub {
            my $frame = shift;
            return $cb->({empty => $frame})
                if $frame->method_frame->isa('Net::AMQP::Protocol::Basic::GetEmpty');
            $self->_push_read_header_and_body('ok', $frame, $cb, $failure_cb);
        },
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub ack {
    my $self = shift;
    my %args = @_;

    return $self if !$self->_check_open(sub {});

    $self->{connection}->_push_write(
        Net::AMQP::Protocol::Basic::Ack->new(
            delivery_tag => 0,
            multiple     => (
                defined $args{delivery_tag} && $args{delivery_tag} != 0 ? 0 : 1
            ),
            %args,
        ),
        $self->{id},
    );

    return $self;
}

sub qos {
    my $self = shift;
    my ($cb, $failure_cb, %args) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Basic::Qos',
        {
            prefetch_count => 1,
            %args,
            prefetch_size  => 0,
            global         => 0,
        },
        'Basic::QosOk', 
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub recover {
    my $self = shift;
    my %args = @_;

    return $self if !$self->_check_open(sub {});

    $self->{connection}->_push_write(
        Net::AMQP::Protocol::Basic::Recover->new(
            requeue => 1,
            %args,
        ),
        $self->{id},
    );

    return $self;
}

sub reject {
    my $self = shift;
    my %args = @_;

    return $self if !$self->_check_open( sub { } );

    $self->{connection}->_push_write(
        Net::AMQP::Protocol::Basic::Reject->new(
            delivery_tag => 0,
            requeue      => 0,
            %args,
        ),
        $self->{id},
    );

    return $self;
}

sub select_tx {
    my $self = shift;
    my ($cb, $failure_cb,) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Tx::Select', {}, 'Tx::SelectOk',
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub commit_tx {
    my $self = shift;
    my ($cb, $failure_cb,) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Tx::Commit', {}, 'Tx::CommitOk',
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub rollback_tx {
    my $self = shift;
    my ($cb, $failure_cb,) = $self->_delete_cbs(@_);

    return $self if !$self->_check_open($failure_cb);

    $self->{connection}->_push_write_and_read(
        'Tx::Rollback', {}, 'Tx::RollbackOk',
        $cb,
        $failure_cb,
        $self->{id},
    );

    return $self;
}

sub push_queue_or_consume {
    my $self = shift;
    my ($frame, $failure_cb,) = @_;

    if ($frame->isa('Net::AMQP::Frame::Method')) {
        my $method_frame = $frame->method_frame;
        if ($method_frame->isa('Net::AMQP::Protocol::Channel::Close')) {
            $self->{connection}->_push_write(
                Net::AMQP::Protocol::Channel::CloseOk->new(),
                $self->{id},
            );
            $self->{_is_open} = 0;
            $self->{_is_active} = 0;
            $self->{_queue}->_flush($frame);
            $self->{_content_queue}->_flush($frame);
            $self->{connection}->delete_channel($self->{id});
            $self->{on_close}->($frame);
            return $self;
        } elsif ($method_frame->isa('Net::AMQP::Protocol::Basic::Deliver')) {
            my $cb = $self->{_consumer_cbs}->{
                $method_frame->consumer_tag
            } || sub {};
            $self->_push_read_header_and_body('deliver', $frame, $cb, $failure_cb);
            return $self;
        } elsif ($method_frame->isa('Net::AMQP::Protocol::Basic::Return')) {
            my $cb = $self->{_return_cbs}->{
                $method_frame->exchange . '_' . $method_frame->routing_key
            } || sub {};
            $self->_push_read_header_and_body('return', $frame, $cb, $failure_cb);
            return $self;
        } elsif ($method_frame->isa('Net::AMQP::Protocol::Channel::Flow')) {
            $self->{_is_active} = $method_frame->active;
            $self->{connection}->_push_write(
                Net::AMQP::Protocol::Channel::FlowOk->new(
                    active => $method_frame->active,
                ),
                $self->{id},
            );
            return $self;
        }
        $self->{_queue}->push($frame);
    } else {
        $self->{_content_queue}->push($frame);
    }

    return $self;
}

sub _push_read_header_and_body {
    my $self = shift;
    my ($type, $frame, $cb, $failure_cb,) = @_;
    my $response = {$type => $frame};
    my $body_size = 0;

    $self->{_content_queue}->get(sub{
        my $frame = shift;

        return $failure_cb->('Received data is not header frame')
            if !$frame->isa('Net::AMQP::Frame::Header');

        my $header_frame = $frame->header_frame;
        return $failure_cb->(
              'Header is not Protocol::Basic::ContentHeader'
            . 'Header was ' . ref $header_frame
        ) if !$header_frame->isa('Net::AMQP::Protocol::Basic::ContentHeader');

        $response->{header} = $header_frame;
        $body_size = $frame->body_size;
    });

    my $body_payload = "";
    my $w_next_frame;
    my $next_frame = sub {
        my $frame = shift;

        return $failure_cb->('Received data is not body frame')
            if !$frame->isa('Net::AMQP::Frame::Body');

        $body_payload .= $frame->payload;

        if (length($body_payload) < $body_size) {
            # More to come
            $self->{_content_queue}->get($w_next_frame);
        }
        else {
            $frame->payload($body_payload);
            $response->{body} = $frame;
            $cb->($response);
        }
    };
    $w_next_frame = $next_frame;
    weaken($w_next_frame);

    $self->{_content_queue}->get($next_frame);

    return $self;
}

sub _delete_cbs {
    my $self = shift;
    my %args = @_;

    my $cb         = delete $args{on_success} || sub {};
    my $failure_cb = delete $args{on_failure} || sub {die @_};

    return $cb, $failure_cb, %args;
}

sub _check_open {
    my $self = shift;
    my ($failure_cb) = @_;

    return 1 if $self->{_is_open};

    $failure_cb->('Channel has already been closed');
    return 0;
}

sub DESTROY {
    my $self = shift;
    $self->close() if defined $self;
    return;
}

1;

1;
__END__

=head1 NAME

AnyEvent::RabbitMQ::Channel - Abstraction of an AMQP channel.

=head1 SYNOPSIS

    my $ch = $rf->open_channel();
    $ch->declare_exchange(exchange => 'test_exchange');

=head1 DESCRIPTION

=head1 METHODS

=head2 declare_exchange (%args)

Declare an exchange (to publish messages to) on the server.

Arguments:

=over

=item on_success

=item on_failure

=item type

Default 'direct'

=item passive

Default 0

=item durable

Default 0

=item auto_delete

Default 0

=item internal

Default 0

=item exchange

The name of the exchange

=back

=head2 delete_exchange

=head2 declare_queue

=head2 bind_queue

Binds a queue to an exchange, with a routing key.

Arguments:

=over

=item queue

The name of the queue to bind

=item exchange

The name of the exchange to bind

=item routing_key

The routing key to bind with

=back

=head2 unbind_queue

=head2 purge_queue

Flushes the contents of a queue.

=head2 delete_queue

Deletes a queue. The queue may not have any active consumers.

=head2 publish

Publish a message to an exchange

Arguments:

=over

=item body

The text body of the message to send.

=item exchange

The name of the exchange to send the message to.

=item routing_key

The routing key with which to publish the message.

=back

=head2 consume

Subscribe to consume messages from a queue.

Arguments:

=over

=item on_consume

Callback called with an argument of the message which has been consumed.

=item consumer_tag

Identifies this consumer, will be auto-generated if you do not provide it, but you must
supply a value if you want to be able to later cancel the subscription.

=item on_success

Callback called if the subscription was successful (before the first message is consumed).

=item on_failure

Callback called if the subscription fails for any reason.

=back

=head2 cancel

Cancel a queue subscription.

Note that the cancellation B<will not> take place at once, and further messages may be
consumed before the subscription is cancelled. No further messages will be
consumed after the on_success callback has been called.

Arguments:

=over

=item consumer_tag

Identifies this consumer, needs to be the value supplied when the queue is initially
consumed from.

=item on_success

Callback called if the subscription was successfully cancelled.

=item on_failure

Callback called if the subscription could not be cancelled for any reason.

=back

=head2 get

Try to get a single message from a queue.

Arguments:

=over

=item queue

Mandatory. Name of the queue to try to receive a message from.

=item on_success

Will be called either with either a message, or, if the queue is empty,
a notification that there was nothing to collect from the queue.

=item on_failure

This callback will be called if an error is signalled on this channel.

=back

=head2 ack

=head2 qos

=head2 recover

=head2 select_tx

=head2 commit_tx

=head2 rollback_tx

=head1 AUTHOR, COPYRIGHT AND LICENSE

See L<AnyEvent::RabbitMQ> for author(s), copyright and license.

=cut


