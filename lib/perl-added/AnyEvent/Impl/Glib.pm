=head1 NAME

AnyEvent::Impl::Glib - AnyEvent adaptor for Glib

=head1 SYNOPSIS

   use AnyEvent;
   use Glib;
  
   # this module gets loaded automatically as required

=head1 DESCRIPTION

This module provides transparent support for AnyEvent. You don't have to
do anything to make Glib work with AnyEvent except by loading Glib before
creating the first AnyEvent watcher.

Glib is probably the most inefficient event loop that has ever seen the
light of the world: Glib not only scans all its watchers (really, ALL of
them, whether I/O-related, timer-related or what not) during each loop
iteration, it also does so multiple times and rebuilds the poll list for
the kernel each time again, dynamically even.

On the positive side, and most importantly, Glib generally works
correctly, no quarrels there.

If you create many watchers (as in: more than two), you might consider one
of the L<Glib::EV>, L<EV::Glib> or L<Glib::Event> modules that map Glib to
other, more efficient, event loops.

This module uses the default Glib main context for all its watchers.

=cut

package AnyEvent::Impl::Glib;

use AnyEvent (); BEGIN { AnyEvent::common_sense }
use Glib ();

our $mainloop = Glib::MainContext->default;

my %io_cond = (
   r => ["in" , "hup"],
   w => ["out", "hup"],
);

sub io {
   my ($class, %arg) = @_;
   
   my $cb = $arg{cb};
   my $fd = fileno $arg{fh};
   defined $fd or $fd = $arg{fh};

   my $source = add_watch Glib::IO
      $fd,
      $io_cond{$arg{poll}},
      sub { &$cb; 1 };

   bless \\$source, $class
}

sub timer {
   my ($class, %arg) = @_;
   
   my $cb   = $arg{cb};
   my $ival = $arg{interval} * 1000;

   my $source; $source = add Glib::Timeout $arg{after} < 0 ? 0 : $arg{after} * 1000,
      $ival ? sub {
                remove Glib::Source $source;
                $source = add Glib::Timeout $ival, sub { &$cb; 1 };
                &$cb;
                0
              }
            : sub { &$cb; 0 };

   bless \\$source, $class
}

sub idle {
   my ($class, %arg) = @_;
   
   my $cb = $arg{cb};
   my $source = add Glib::Idle sub { &$cb; 1 };
   bless \\$source, $class
}

sub DESTROY {
   remove Glib::Source $${$_[0]};
}

sub one_event {
   $mainloop->iteration (1);
}

sub loop {
   # hackish, but we do not have a mainloop, just a maincontext
   $mainloop->iteration (1) while 1;
}

1;

=head1 SEE ALSO

L<AnyEvent>, L<Glib>.

=head1 AUTHOR

 Marc Lehmann <schmorp@schmorp.de>
 http://home.schmorp.de/

=cut

