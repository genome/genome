
package Genome::Model::Command::Services::WebApp::Starman::Server;

use base qw( Starman::Server );
use Socket qw(IPPROTO_TCP TCP_NODELAY SOL_SOCKET SO_REUSEADDR);

sub process_request {
    my $self = shift;
    my $conn = $self->{server}->{client};

    if ($conn->NS_proto eq 'TCP') {
        setsockopt($conn, IPPROTO_TCP, TCP_NODELAY, 1)
            or die $!;
        setsockopt($conn, SOL_SOCKET, SO_REUSEADDR, 1)
            or die $!;
    }

    while ( $self->{client}->{keepalive} ) {
        last if !$conn->connected;

        # Read until we see all headers
        last if !$self->_read_headers;

        my $env = {
            REMOTE_ADDR     => $self->{server}->{peeraddr},
            REMOTE_HOST     => $self->{server}->{peerhost} || $self->{server}->{peeraddr},
            SERVER_NAME     => $self->{server}->{sockaddr}, # XXX: needs to be resolved?
            SERVER_PORT     => $self->{server}->{sockport},
            SCRIPT_NAME     => '',
            'psgi.version'      => [ 1, 1 ],
            'psgi.errors'       => *STDERR,
            'psgi.url_scheme'   => 'http',
            'psgi.nonblocking'  => Plack::Util::FALSE,
            'psgi.streaming'    => Plack::Util::TRUE,
            'psgi.run_once'     => Plack::Util::FALSE,
            'psgi.multithread'  => Plack::Util::FALSE,
            'psgi.multiprocess' => Plack::Util::TRUE,
            'psgix.io'          => $conn,
            'psgix.input.buffered' => Plack::Util::TRUE,
        };

        # Parse headers
        my $reqlen = parse_http_request(delete $self->{client}->{headerbuf}, $env);
        if ( $reqlen == -1 ) {
            # Bad request
            DEBUG && warn "[$$] Bad request\n";
            $self->_http_error(400, { SERVER_PROTOCOL => "HTTP/1.0" });
            last;
        }

        # Initialize PSGI environment
        # Determine whether we will keep the connection open after the request
        my $connection = delete $env->{HTTP_CONNECTION};
        my $proto = $env->{SERVER_PROTOCOL};
        if ( $proto && $proto eq 'HTTP/1.0' ) {
            if ( $connection && $connection =~ /^keep-alive$/i ) {
                # Keep-alive only with explicit header in HTTP/1.0
                $self->{client}->{keepalive} = 1;
            }
            else {
                $self->{client}->{keepalive} = 0;
            }
        }
        elsif ( $proto && $proto eq 'HTTP/1.1' ) {
            if ( $connection && $connection =~ /^close$/i ) {
                $self->{client}->{keepalive} = 0;
            }
            else {
                # Keep-alive assumed in HTTP/1.1
                $self->{client}->{keepalive} = 1;
            }

            # Do we need to send 100 Continue?
            if ( $env->{HTTP_EXPECT} ) {
                if ( $env->{HTTP_EXPECT} eq '100-continue' ) {
                    syswrite $conn, 'HTTP/1.1 100 Continue' . $CRLF . $CRLF;
                    DEBUG && warn "[$$] Sent 100 Continue response\n";
                }
                else {
                    DEBUG && warn "[$$] Invalid Expect header, returning 417\n";
                    $self->_http_error( 417, $env );
                    last;
                }
            }

            unless ($env->{HTTP_HOST}) {
                # No host, bad request
                DEBUG && warn "[$$] Bad request, HTTP/1.1 without Host header\n";
                $self->_http_error( 400, $env );
                last;
            }
        }

        unless ($self->{options}->{keepalive}) {
            DEBUG && warn "[$$] keep-alive is disabled. Closing the connection after this request\n";
            $self->{client}->{keepalive} = 0;
        }

        $self->_prepare_env($env);

        # Run PSGI apps
        my $res = Plack::Util::run_app($self->{app}, $env);

        if (ref $res eq 'CODE') {
            $res->(sub { $self->_finalize_response($env, $_[0]) });
        } else {
            $self->_finalize_response($env, $res);
        }

        DEBUG && warn "[$$] Request done\n";

        if ( $self->{client}->{keepalive} ) {
            # If we still have data in the input buffer it may be a pipelined request
            if ( $self->{client}->{inputbuf} ) {
                if ( $self->{client}->{inputbuf} =~ /^(?:GET|HEAD)/ ) {
                    if ( DEBUG ) {
                        warn "Pipelined GET/HEAD request in input buffer: " 
                            . dump( $self->{client}->{inputbuf} ) . "\n";
                    }

                    # Continue processing the input buffer
                    next;
                }
                else {
                    # Input buffer just has junk, clear it
                    if ( DEBUG ) {
                        warn "Clearing junk from input buffer: "
                            . dump( $self->{client}->{inputbuf} ) . "\n";
                    }

                    $self->{client}->{inputbuf} = '';
                }
            }

            DEBUG && warn "[$$] Waiting on previous connection for keep-alive request...\n";

            my $sel = IO::Select->new($conn);
            last unless $sel->can_read(1);
        }
    }

    DEBUG && warn "[$$] Closing connection\n";

    if ($self->{options}->{single_request}) {
        DEBUG && warn "[$$] Exiting\n";
        exit;
    }
}

1;
