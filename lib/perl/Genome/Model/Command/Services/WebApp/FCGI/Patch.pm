package Genome::Model::Command::Services::WebApp::FCGI::Patch;

use strict;
use Socket;
#use File::FDpasser;
use FCGI;
use POSIX;
use Time::HiRes qw(usleep);

#print "loaded fcgifix\n";

#
# This module monkeypatches over the FCGI::OpenSocket function
# so we can connect to the old process and tell it we're ready
# to take over the listener
#
#
{
    no strict 'refs';

    *{"FCGI::__OpenSocket"} = \&FCGI::OpenSocket;
}

our $FD;

$SIG{USR1} = sub {
    POSIX::close($FD);
};

sub open ($$) {
    if ( exists $ENV{RESTARTSOCK} ) {
        my $s = $ENV{RESTARTSOCK};

        #warn "connecting to $s to get listener: @_\n";

        socket( PEER, PF_UNIX, SOCK_STREAM, 0 ) or die "socket: $!";
        my $retry = 0;
        while ( $retry++ < 3 ) {
            connect( PEER, sockaddr_un($s) ) and last;
            sleep 1;
        }

        my $r = <PEER>;
        close(PEER);

        usleep(2e5);
    }
    $FD = FCGI::__OpenSocket(@_);

    #warn "got fileno: $FD\n";
    return $FD;
}

*FCGI::OpenSocket = \&open;

1;
