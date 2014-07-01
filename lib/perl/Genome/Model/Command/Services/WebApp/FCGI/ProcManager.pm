
package Genome::Model::Command::Services::WebApp::FCGI::ProcManager;

use strict;
use warnings;
use base qw(FCGI::ProcManager);

# I spent a day trying to make this work with just basic system calls
# like select() and signal handlers.  It was so hard to understand
# without reimplementing pm_manage()  that I gave up and am now using
# AnyEvent.
#
use EV;
use AnyEvent;
use AnyEvent::Socket;
use AnyEvent::Handle;
use AnyEvent::Util;
use POSIX qw(setsid);
#use File::FDpasser;

sub _disabled_pm_pre_dispatch {
    my $self = shift;

#    $self->pm_notify('pre_dispatch');

    $self->SUPER::pm_pre_dispatch(@_);
}

sub pm_post_dispatch {
    my $self = shift;

#    $self->pm_notify('post_dispatch');

    exit;
#    UR::Context->rollback;
#    UR::Context->clear_cache;

    $self->SUPER::pm_post_dispatch(@_);
}

my $sock = "/tmp/gvsock";
my %w = ();
my $restart_cv;
my $safe_cleanup = 0;
sub pm_wait {
    my $self = shift;

    # see the AnyEvent docs for the AE api, its faster with most event loops

    foreach my $kp (keys %{ $self->{PIDS}}) {
        next if exists($w{$kp});
        $w{$kp} = AE::child($kp,sub {
            my $pid;
            ($pid, $?) = @_;

            delete $w{$pid};
            delete $self->{PIDS}->{$pid} and
                $self->pm_notify("server (pid $pid) exited with status $?");

            EV::unloop();
        });
    }

    # this is assuming imp-apipe
    my $restart_command = '/usr/local/sbin/genome_view';

    ### wait for HUP then start the command and a listener.

    if (!exists $w{'HUP'}) {
        $w{'HUP'} = AE::signal('HUP', sub {
            $SIG{'HUP'} = 'IGNORE';
            warn "Restart triggered\n";

            $restart_cv = run_cmd(
                $restart_command,
                '>' => '/dev/null',
                '2>' => '/dev/null',
                '<' => '/dev/null',
                close_all => 1,
                on_prepare => sub {
                    $ENV{'RESTARTSOCK'} = $sock;
                    chdir '/';
                    setsid;
            });
            $restart_cv->cb(sub {
                warn "restart script unexpected exit: " . $_[0]->recv;
                undef $restart_cv;
                delete $w{'listener'};
                delete $w{'wtr'};
                unlink $sock;
                delete $w{'HUP'};

                if ($self->pid_fname) {
                    open(P, ">" . $self->pid_fname) or return;
                    print P $$ . "\n";
                    close P;
                }
            });

            $w{'listener'} = tcp_server "unix/", $sock, sub {
                my ($fh) = @_;

                delete $w{'listener'}; ## only accept the first connection
                unlink $sock;

                $w{'wtr'} = AE::io($fh, 1, sub {
                    $fh->print('1');
#                    send_file(fileno($fh),$FCGIfix::FD);
                    close($fh);

                    $safe_cleanup = 1;
                    EV::unloop;
                });

            };

            $w{'HUP'} = 0;
        });
    }

    EV::loop();
    if ($safe_cleanup) {
        $SIG{USR1}->();
        kill 'USR1', keys %{$self->{PIDS}} if $self->{PIDS};

        while (keys %{$self->{PIDS}}) {
            $self->SUPER::pm_wait();
        }

        $self->pm_notify("server all children dead, exiting");

        exit 0;
    }
}

1;
