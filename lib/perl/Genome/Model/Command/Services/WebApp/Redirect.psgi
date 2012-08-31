use strict;
use warnings;

use Data::Dumper;

sub {
    my ( $env, $relative, $code ) = @_;
    $code ||= 302;

    my $static;
    if ( $relative =~ /^https?:/ ) {
        $static = $relative;
    } else {
        my $http_host =
          defined $env->{'HTTP_HOST'}
          ? $env->{'HTTP_HOST'}
          : do {
            my $server_name = $env->{'SERVER_NAME'};
            my $server_port = $env->{'SERVER_PORT'};
            $server_name . ( $server_port != 80 ? ":" . $server_port : '' );
          };

        $static = $env->{'psgi.url_scheme'} . '://' . $http_host . $relative;
    }

    return [ $code, [ 'Location', $static ], [] ];
};
