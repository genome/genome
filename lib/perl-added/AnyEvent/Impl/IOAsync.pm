=head1 NAME

AnyEvent::Impl::IOAsync - AnyEvent adaptor for IO::Async

=head1 SYNOPSIS

  use AnyEvent;

  use IO::ASync::Loop;
  use AnyEvent::Impl::IOAsync;

  my $loop = new IO::Async::Loop;

  AnyEvent::Impl::IOAsync::set_loop $loop;

=head1 DESCRIPTION

This module provides support for IO::Async as AnyEvent backend. Due to
the rather sad state of L<IO::Async>, support is only available partially
(only timers and I/O watchers are supported, signals and child watchers
are emulated by AnyEvent itself (fighting with IO::Async, so you cannot
use both), idle watchers are being emulated, I/O watchers need to dup
their fh.

=head1 PROBLEMS WITH IO::Async

There have been multiple attempts at providing an AnyEvent interface to
IO::Async, and the effort is ongoing. Supporting IO::Async is hard. Here's
why:

=over 4

=item IO::Async integration cannot be automatic - no default loop

IO::Async doesn't offer an interface suitable for independent usage of
event sources: there is no standard way to share the main event loop
between modules - modules have to somehow agree on how to do this.

For AnyEvent to work with IO::Async, the IO::Async main program has to
call C<AnyEvent::Impl::IOAsync::set_loop> with the C<IO::Async::Loop>
object that AnyEvent is to use, see the SYNOPSIS section for an example.

It is possible to get a copy of the loop by reading
C<$AnyEvent::Impl::IOAsync::LOOP>, so AnyEvent could be used as a central
place to store the default L<IO::Async::Loop> object, also for other
modules, if there is any need for that.

=item Broken child watchers

IO::Async is the only module that requires you to call a special function
before even forking your child program, while the AnyEvent API itself
(which works with other event loops) works as long as it has been
initialised (AnyEvent's own pure perl implementation doesn't even suffer
from these limitations and just works).

Worse, IO::Async does not let you install multiple child watchers, does
not let you watch for any child, and apparently an interface to unregister
child watchers has been forgotten as well.

As a result, AnyEvent::Impl::IOAsync has to fall back on it's own child
management, which makes it impossible to watch for child processes via
both AnyEvent and IO::Async. Hooking and Patching IO::Async has been
considered, but is considerable work.

=item No support for multiple watchers per event

In most (all? documentation?) cases you cannot have multiple watchers
for the same event (what's the point of having all these fancy notifier
classes when you cannot have multiple notifiers for the same event? That's
like only allowing one timer per second or so...).

This makes signal watchers almost useless (You could just hook them
yourself, you can't share any of them, as would make sense for
e.g. SIGTERM, SIGTSTP, SIGPWR, SIGUSR1 etc.).

As a result, AnyEvent falls back to it's own signal handling (it is
pointless to somehow share the IO::Async watcher, as it doesn't matter if
AnyEvent blocks the signal via IO::Async or directly, and AnyEvents signal
handling is race-free).

For I/O watchers, AnyEvent has to dup() every file handle, as IO::Async
fails to support the same or different file handles pointing to the same
fd (this is at least documented, but why not fix it instead?).

=back

Apart from these fatal flaws, there are a number of unpleasent properties
that just need some mentioning:

=over 4

=item Confusing and misleading name

Another rather negative point about this module family is its name,
which is deeply confusing: Despite the "async" in the name, L<IO::Async>
only does I<synchronous> I/O, there is nothing "asynchronous" about it
whatsoever (when I first heard about it, I thought, "wow, a second async
I/O module, what does it do compared to L<IO::AIO>", and was somehow set
back when I learned that the only "async" aspect of it is the name).

=item Inconsistent, incomplete and convoluted API

Implementing AnyEvent's rather simple timers on top of IO::Async's timers
was a nightmare (try implementing a timer with configurable interval and
delay value...).

How to actually get I/O events in L<IO::Async::Handle> is undocumented:
read events are apparently automatic, for write events, you have to
explicitly request C<want_writeready>, and specifying C<want_readready> is
apparently a usage bug (it doesn't exist). All this must be deduced from
reading the sources.

You can't stop child watchers. Even reading the sources I found no way to
stop them. It must have been forgotten.

The method naming is chaotic: C<watch_child> creates a child watcher,
but C<watch_io> is an internal method; C<detach_signal> removes a signal
watcher, but C<detach_child> forks a subprocess and so on).

IO::Async has weird checks - passing in a callable reference is sometimes
forbidden (of course, this is checked on every invocation, not when the
callback is registered, so you have no idea where in your code you passed
it in), as the code checks explicitly for code references, disallowing
callable objects.

=item Unpleasant surprises on GNU/Linux

When you develop your program on FreeBSD and run it on GNU/Linux, you
might have unpleasant surprises, as IO::Async::Loop will by default use
L<IO::Async::Loop::Epoll>, which is incompatible with C<fork>, so your
network server will run into spurious and very hard to debug problems
under heavy load, as IO::Async forks a lot of processes, e.g. for DNS
resolution. It would be better if IO::Async would only load "safe"
backends by default (or fix the epoll backend to work in the presence of
fork, which admittedly is hard - EV does it for you, and also does not use
unsafe backends by default).

=item Exiting considered harmful

   (in cleanup) Can't call method "parent" on an undefined value
      at IO/Async/Loop.pm line 297 during global destruction.

IO::Async just hates global destruction. Calling C<exit> will easily give
you one such line per watcher.

The problem is that L<IO::Async::Loop> is itself not warning-free, but
actually enables warnings for itself.

(Ok, the real bug is of course perl's broken mark & sweep garbage
collector that corrupts data structures).

=back

On the positive side, performance with IO::Async is quite good even in my
very demanding eyes.

=cut

package AnyEvent::Impl::IOAsync;

use AnyEvent (); BEGIN { AnyEvent::common_sense }

use Time::HiRes;
use Scalar::Util;

BEGIN {
   # IO::Async enables warnings but is itself not warning-free, we
   # try our best to silence it.
   require warnings;
   local *warnings::import = sub { };
   require IO::Async::Loop;
}

use IO::Async::Handle;

our $VERSION = $AnyEvent::VERSION;
our $LOOP;

sub set_loop($) {
   $LOOP = $_[0];
   #$LOOP->enable_childmanager;
}

sub timer {
   my ($class, %arg) = @_;
   
   # IO::Async has problems with overloaded objects
   my $cb = $arg{cb};

   my $id;

   if (my $ival = $arg{interval}) {
      my $ival_cb; $ival_cb = sub {
         $id = $LOOP->enqueue_timer (delay => $ival, code => $ival_cb);
         &$cb;
      };
      $id = $LOOP->enqueue_timer (delay => $arg{after}, code => $ival_cb);

      # we have to weaken afterwards, but when enqueue dies, we have a memleak.
      # still, we do anything for speed...
      Scalar::Util::weaken $ival_cb;

   } else {
      $id = $LOOP->enqueue_timer (delay => $arg{after}, code => sub { &$cb });
   }

   bless \\$id, "AnyEvent::Impl::IOAsync::timer"
}

sub AnyEvent::Impl::IOAsync::timer::DESTROY {
   $LOOP->cancel_timer (${${$_[0]}});
}

sub io {
   my ($class, %arg) = @_;

   # we need to dup(), as IO::Async only allows one watcher per
   # underlying fd.
   my ($fh) = AnyEvent::_dupfh $arg{poll}, $arg{fh};

   # there is no want_readready, but want_writeready is obligatory :/
   # I just love undocumented and illogical interfaces...
   my $id = new IO::Async::Handle
      $arg{poll} eq "r"
         ? (read_handle  => $fh, on_read_ready  => $arg{cb})
         : (write_handle => $fh, on_write_ready => $arg{cb}, want_writeready => 1),
   ;
   $LOOP->add ($id);

   bless \\$id, "AnyEvent::Impl::IOAsync::io"
}

sub AnyEvent::Impl::IOAsync::io::DESTROY {
   $LOOP->remove (${${$_[0]}});
}

#sub signal {
#   my ($class, %arg) = @_;
#
#   my $signal = $arg{signal};
#
#   my $id = $LOOP->attach_signal ($arg{signal}, $arg{cb});
#   bless [$signal, $id], "AnyEvent::Impl::IOAsync::signal";
#}
#
#sub AnyEvent::Impl::IOAsync::signal::DESTROY {
#   $LOOP->detach_signal (@{ $_[0] });
#}

sub one_event {
   $LOOP->loop_once;
}

sub loop {
   $LOOP->loop_forever;
}

1;

=head1 SEE ALSO

L<AnyEvent>, L<EV>.

=head1 AUTHOR

 Marc Lehmann <schmorp@schmorp.de>
 http://home.schmorp.de/

=cut

