package AnyEvent;
sub CYGWIN () { 0 }
sub WIN32 () { 0 }
sub F_SETFL () { 4 }
sub F_SETFD () { 2 }
sub O_NONBLOCK () { 2048 }
sub FD_CLOEXEC () { 1 }
package AnyEvent::Util;
sub WSAEINVAL () { -1e+99 }
sub WSAEWOULDBLOCK () { -1e+99 }
sub WSAEINPROGRESS () { -1e+99 }
sub _AF_INET6 () { 10 }
1;
