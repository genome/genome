=head1 NAME

AnyEvent::SNMP - adaptor to integrate Net::SNMP into AnyEvent.

=head1 SYNOPSIS

 use AnyEvent::SNMP;
 use Net::SNMP;

 # just use Net::SNMP and AnyEvent as you like:

 # use a condvar to transfer results, this is
 # just an example, you can use a naked callback as well.
 my $cv = AnyEvent->condvar;

 # ... start non-blocking snmp request(s)...
 Net::SNMP->session (-hostname => "127.0.0.1",
                     -community => "public",
                     -nonblocking => 1)
          ->get_request (-callback => sub { $cv->send (@_) });

 # ... do something else until the result is required
 my @result = $cv->wait;

=head1 DESCRIPTION

This module implements an alternative "event dispatcher" for Net::SNMP,
using AnyEvent as a backend.

This integrates Net::SNMP into AnyEvent: You can make non-blocking
Net::SNMP calls and as long as other parts of your program also use
AnyEvent (or some event loop supported by AnyEvent), they will run in
parallel.

Also, the Net::SNMP scheduler is very inefficient with respect to both CPU
and memory usage. Most AnyEvent backends (including the pure-perl backend)
fare much better than the Net::SNMP dispatcher.

A potential disadvantage is that replacing the dispatcher is not at all
a documented thing to do, so future changes in Net::SNP might break this
module (or the many similar ones).

This module does not export anything and does not require you to do
anything special apart from loading it I<before doing any non-blocking
requests with Net::SNMP>. It is recommended but not required to load this
module before C<Net::SNMP>.

=head1 GLOBAL VARIABLES

=over 4

=item $AnyEvent::SNMP::MAX_OUTSTANDING (default: C<50>, dynamic)

=item AnyEvent::SNMP::set_max_outstanding $new_value

Use this package variable to restrict the number of outstanding SNMP
requests at any point in time.

Net::SNMP is very fast at creating and sending SNMP requests, but much
slower at parsing (big, bulk) responses. This makes it easy to request a
lot of data that can take many seconds to parse.

In the best case, this can lead to unnecessary delays (and even time-outs,
as the data has been received but not yet processed) and in the worst
case, this can lead to packet loss, when the receive queue overflows and
the kernel can no longer accept new packets.

To avoid this, you can (and should) limit the number of outstanding
requests to a number low enough so that parsing time doesn't introduce
noticable delays.

Unfortunately, this number depends not only on processing speed and load
of the machine running Net::SNMP, but also on the network latency and the
speed of your SNMP agents.

AnyEvent::SNMP tries to dynamically adjust this number dynamically upwards
and downwards.

Increasing C<$MAX_OUTSTANDING> will not automatically use the
extra request slots. To increase C<$MAX_OUTSTANDING> and make
C<AnyEvent::SNMP> make use of the extra paralellity, call
C<AnyEvent::SNMP::set_max_outstanding> with the new value, e.g.:

   AnyEvent::SNMP::set_max_outstanding 500;

Although due to the dynamic adjustment, this might have little lasting
effect.

Note that you can use L<Net::SNMP::XS> to speed up parsing of responses
considerably.

=item $AnyEvent::SNMP::MIN_RECVQUEUE (default: C<8>)

=item $AnyEvent::SNMP::MAX_RECVQUEUE (default: C<64>)

These values specify the minimum and maximum receive queue length (in
units of one response packet).

When AnyEvent::SNMP handles $MAX_RECVQUEUE or more packets per iteration
it will reduce $MAX_OUTSTANDING. If it handles less than $MIN_RECVQUEUE,
it increases $MAX_OUTSTANDING.

This has the result of adjusting the number of outstanding requests so that
the recv queue is between the minimum and maximu, usually.

This algorithm works reasonably well as long as the responses, response
latencies and processing times are the same size per packet on average.

=back

=head1 COMPATIBILITY

This module may be used as a drop in replacement for the
Net::SNMP::Dispatcher in existing programs. You can still call
C<snmp_dispatcher> to start the event-loop, but then you loose the benefit
of mixing Net::SNMP events with other events.

   use AnyEvent::SNMP;
   use Net::SNMP;

   # just use Net::SNMP as before

   # ... start non-blocking snmp request(s)...
   Net::SNMP->session (
         -hostname    => "127.0.0.1",
         -community   => "public",
         -nonblocking => 1,
      )->get_request (-callback => sub { ... });

   snmp_dispatcher;

=cut

package AnyEvent::SNMP;

no warnings;
use strict qw(subs vars);

# it is possible to do this without loading
# Net::SNMP::Dispatcher, but much more awkward.
use Net::SNMP::Dispatcher;

sub Net::SNMP::Dispatcher::instance {
   AnyEvent::SNMP::
}

use Net::SNMP ();
use AnyEvent ();

our $VERSION = '1.0';

$Net::SNMP::DISPATCHER = instance Net::SNMP::Dispatcher;

our $MESSAGE_PROCESSING = $Net::SNMP::Dispatcher::MESSAGE_PROCESSING;

our $BUSY;
our $DONE; # finished all jobs
our @TRANSPORT; # fileno => [count, watcher]
our @QUEUE;
our $MAX_OUTSTANDING = 50;
our $MIN_RECVQUEUE   =  8;
our $MAX_RECVQUEUE   = 64;

sub kick_job;

sub _send_pdu {
   my ($pdu, $retries) = @_;

   # mostly copied from Net::SNMP::Dispatch

   # Pass the PDU to Message Processing so that it can
   # create the new outgoing message.
   my $msg = $MESSAGE_PROCESSING->prepare_outgoing_msg ($pdu);

   if (!defined $msg) {
      --$BUSY;
      kick_job;
      # Inform the command generator about the Message Processing error.
      $pdu->status_information ($MESSAGE_PROCESSING->error);
      return; 
   }

   # Actually send the message.
   if (!defined $msg->send) {
      $MESSAGE_PROCESSING->msg_handle_delete ($pdu->msg_id)
         if $pdu->expect_response;

      # A crude attempt to recover from temporary failures.
      if ($retries-- > 0 && ($!{EAGAIN} || $!{EWOULDBLOCK} || $!{ENOSPC})) {
         my $retry_w; $retry_w = AE::timer $pdu->timeout, 0, sub {
            undef $retry_w;
            _send_pdu ($pdu, $retries);
         };
      } else {
         --$BUSY;
         kick_job;
      }

      # Inform the command generator about the send() error.
      $pdu->status_information ($msg->error);
      return;
   }

   # Schedule the timeout handler if the message expects a response.
   if ($pdu->expect_response) {
      my $transport = $msg->transport;
      my $fileno    = $transport->fileno;

      # register the transport
      unless ($TRANSPORT[$fileno][0]++) {
         $TRANSPORT[$fileno][1] = AE::io $transport->socket, 0, sub {
            for my $count (1..$MAX_RECVQUEUE) { # handle up to this many requests in one go
               # Create a new Message object to receive the response
               my ($msg, $error) = Net::SNMP::Message->new (-transport => $transport);

               if (!defined $msg) {
                  die sprintf 'Failed to create Message object [%s]', $error;
               }

               # Read the message from the Transport Layer
               if (!defined $msg->recv) {
                  if ($transport->connectionless) {
                     # if we handled very few replies and we have queued work, try
                     # to increase the parallelity as we probably can handle more.
                     if ($count < $MIN_RECVQUEUE && @QUEUE) {
                        ++$MAX_OUTSTANDING;
                        kick_job;
                     }
                  } else {
                     # for some reason, connected-oriented transports seem to need this
                     delete $TRANSPORT[$fileno]
                        unless --$TRANSPORT[$fileno][0];
                  }

                  $msg->error;
                  return;
               }

               # For connection-oriented Transport Domains, it is possible to
               # "recv" an empty buffer if reassembly is required.
               if (!$msg->length) {
                  return;
               }

               # Hand the message over to Message Processing.
               if (!defined $MESSAGE_PROCESSING->prepare_data_elements ($msg)) {
                  $MESSAGE_PROCESSING->error;
                  return;
               }

               # Set the error if applicable.
               $msg->error ($MESSAGE_PROCESSING->error) if $MESSAGE_PROCESSING->error;

               # Notify the command generator to process the response.
               $msg->process_response_pdu; 

               # Cancel the timeout.
               my $rtimeout_w = $msg->timeout_id;
               if ($$rtimeout_w) {
                  undef $$rtimeout_w;

                  --$BUSY;
                  kick_job;

                  unless (--$TRANSPORT[$fileno][0]) {
                     delete $TRANSPORT[$fileno];
                     return;
                  }
               }
            }

            # when we end up here, we successfully handled $MAX_RECVQUEUE
            # replies in one iteration, so assume we are overloaded
            # and reduce the amount of parallelity.
            $MAX_OUTSTANDING = (int $MAX_OUTSTANDING * 0.95) || 1;
         };
      }

      $msg->timeout_id (\(my $rtimeout_w =
         AE::timer $pdu->timeout, 0, sub {
            my $rtimeout_w = $msg->timeout_id;
            if ($$rtimeout_w) {
               undef $$rtimeout_w;
               delete $TRANSPORT[$fileno]
                  unless --$TRANSPORT[$fileno][0];
            }

            if ($retries--) {
               _send_pdu ($pdu, $retries);
            } else {
               $MESSAGE_PROCESSING->msg_handle_delete ($pdu->msg_id);
               $pdu->status_information ("No response from remote host '%s'", $pdu->hostname);

               --$BUSY;
               kick_job;
            }
         })
      ); 
   } else {
     --$BUSY;
     kick_job;
   }
}

sub kick_job {
   while ($BUSY < $MAX_OUTSTANDING) {
      my $pdu = shift @QUEUE
         or last;

      ++$BUSY;
      _send_pdu $pdu, $pdu->retries;
   }

   $DONE and $DONE->() unless $BUSY;
}

sub send_pdu($$$) {
   my (undef, $pdu, $delay) = @_;

   # $delay is not very sensibly implemented by AnyEvent::SNMP,
   # but apparently it is not a very sensible feature.
   if ($delay > 0) {
      ++$BUSY;
      my $delay_w; $delay_w = AE::timer $delay, 0, sub {
         undef $delay_w;
         push @QUEUE, $pdu;
         --$BUSY;
         kick_job;
      };
      return 1;
   }

   push @QUEUE, $pdu;
   kick_job;

   1
}

sub activate($) {
   while ($BUSY) {
      $DONE = AE::cv;
      $DONE->recv;
      undef $DONE;
   }
}

sub one_event($) {
   # should not ever be used
   AnyEvent->one_event; #d# todo
}

sub set_max_outstanding($) {
   $MAX_OUTSTANDING = $_[0];
   kick_job;
}

=head1 SEE ALSO

L<AnyEvent>, L<Net::SNMP>, L<Net::SNMP::XS>, L<Net::SNMP::EV>.

=head1 AUTHOR

 Marc Lehmann <schmorp@schmorp.de>
 http://home.schmorp.de/

=cut

1

