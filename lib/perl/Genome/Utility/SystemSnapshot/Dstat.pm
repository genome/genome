package Genome::Utility::SystemSnapshot::Dstat;

use Storable;
use Data::Dumper;

use base qw(Genome::Utility::SystemSnapshot);
use Genome::Utility::AsyncFileSystem;

use strict;
use warnings;

sub _sleep {
  # Use an event timer to sleep within the event loop
  my $self = shift;
  my $seconds = shift;
  my $sleep_cv = AnyEvent->condvar;
  my $w = AnyEvent->timer(
      after => $seconds,
      cb => $sleep_cv
  );
  $sleep_cv->recv;
}

sub pre_run {
  my $self = shift;
  my $cmd = "/gsc/var/gsc/systems/blades/dstat/dstat";
  if (! -x $cmd) {
    warn "comand not found: $cmd: NOT profiling";
    return;
  }
  #my $args = "-Tcldigmnpsy --ipc --lock --tcp --udp --unix --output $self->{metrics}";
  my $args;
  if ($ENV{DSTAT_ARGS}) {
      $args = "--output $self->{metrics} ";
      $args .= $ENV{DSTAT_ARGS};
  } else {
      $args = "-Tcldigmnpsy --ipc --lock --tcp --udp --unix --output $self->{metrics}";
  }

  # Add a command to a condition variable event loop.
  # This gets started by our caller's use of $cmd_cv->recv;
  $self->{cv} = Genome::Utility::AsyncFileSystem->shellcmd(
      cmd => "$cmd $args",
      '>' => '/dev/null',
      '2>' => '/dev/null',
      '$$' => \( $self->{pid} ),
      close_all => 1,
      allow_failed_exit_code => 1,
  );
  # Give cmd time to start up
  $self->_sleep(2);
}

sub run {
  my $self = shift;
  my $cmd = shift;
  my $args = shift;

  my $cmd_cv = Genome::Utility::AsyncFileSystem->shellcmd(
    '>' => $self->{output},
    '2>' => $self->{errors},
    cmd => "$cmd $args"
  );

  # This begins the event loop that runs both the snapshotter and the cmd
  $cmd_cv->recv;
}

sub post_run {
  my $self = shift;

  # Now that cmd_cv->cmd status is true, we're back, and we send SIGTERM to collectl's pid.
  kill 15, $self->{pid};
  # Now recv on that condition variable, which will catch the signal and exit.
  # We wrap in eval and examine $@ to ensure we catch the signal we sent, but we can still
  # observe any unexpected events.
  eval {
    $self->{cv}->recv;
  };
  if ($@ && $@ !~ /^COMMAND KILLED\. Signal 15/) {
    # unexpected death message from shellcmd.
    # Don't die here, we don't want shutting down the profiler to interrupt
    # the job being profiled.  Just warn.
    warn "unexpected death of profiler: $@";
  }
}

sub report {
  my $self = shift;
  my $metrics = {};
  # dstat doesn't report key/value metrics, return nothing.
  return $metrics;
}

1;
