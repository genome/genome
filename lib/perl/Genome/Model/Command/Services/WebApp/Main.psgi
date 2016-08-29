#!/usr/bin/env genome-perl

use Genome;
use Web::Simple 'Genome::Model::Command::Services::WebApp::Main';

package Genome::Model::Command::Services::WebApp::Main;

use strict;

use File::Basename;
use Plack::Util;
use Plack::Builder;

our $always_memcache = Genome::Config::get('view_cache');

our $psgi_path;
eval {
    ## Genome is probably not loaded, thats okay, we can just use __FILE__
    $psgi_path = Genome::Model::Command::Services::WebApp->psgi_path;
};
unless (defined $psgi_path) {
    $psgi_path = File::Basename::dirname(__FILE__);
}

my @psgi = qw(
    Rest.psgi
    Redirect.psgi
    404Handler.psgi
    Dump.psgi
    Cache.psgi
    Info.psgi
    File.psgi
);
our %app = map { $_ => load_app($_) } @psgi;

## Utility functions
sub load_app {
    Plack::Util::load_psgi( $psgi_path . '/' . shift );
}

sub redispatch_psgi {
    my ( $psgi_app, $env, @args ) = @_;

    if (ref($env) eq 'HASH' && $env->{'REMOTE_USER'} ) {
        $ENV{'REMOTE_USER'} = $env->{'REMOTE_USER'};
    } else {
        # dev server doesnt have REMOTE_USER set
        $ENV{'REMOTE_USER'} = Genome::Sys->username();
    }

    my $response = $psgi_app->( $env, @args );
    return sub  { $response };
}

sub redirect_to {
    my ($path, $request) = @_;
    $request->{REDIRECT_URI} = $path;
    redispatch_psgi( $app{'Redirect.psgi'}, $request);
}

## Web::Simple dispatcher for all apps
sub dispatch_request {
    ## make 404's pretty by sending them to 404Handler.psgi

    my ($self, $request) = @_;

    if (defined $request->{SCRIPT_NAME} && $request->{SCRIPT_NAME} ne ''  && ($request->{PATH_INFO} ne "/" && $request->{PATH_INFO} ne ""))  {
        $request->{PATH_INFO} = $request->{SCRIPT_NAME} . $request->{PATH_INFO};
    }
    if ($request->{PATH_INFO} eq "") {
        $request->{PATH_INFO} = "/";
    }

    sub (GET) {
        response_filter {
            my $resp = $_[0];
            if ( ref($resp) eq 'ARRAY' && $resp->[0] == 404 ) {
                return redispatch_psgi( $app{'404Handler.psgi'}, $resp->[2] );
            } elsif ( ref($resp) eq 'ARRAY' && $resp->[0] == 500 ) {
                return redispatch_psgi( $app{'404Handler.psgi'}, $resp->[2] );
            }
            return $resp;
        };

    },
    sub (/view/whoami) {
        my $u = $request->{'REMOTE_USER'} || Genome::Sys->username();
        my $sys_user = Genome::Sys::User->get(username => $u);
        my $c = { username => $u, id => $sys_user->id };
        use JSON;
        my $body = to_json( $c, { ascii => 1, allow_nonref => 1, });
        return [ 200, [ 'Content-type'   => "text/plain" ], [$body] ];
    },
    sub (/view/x/...) {
        redispatch_psgi($app{'File.psgi'}, $request);
    },
    sub (/view/debug) {
        redispatch_psgi($app{'Info.psgi'}, $request);
    },
    sub (/res/** + .*) {
        my $new_path = "/genome/resource.html/$_[1].$_[2]";
        my %new_params = ( %{$_[3]}, PATH_INFO => $new_path, REQUEST_URI => $new_path );
        redispatch_psgi($app{'Rest.psgi'}, \%new_params);
    },
      ## send /view without a trailing slash to /view/
      ## although thats probably a 404
      sub (/view) {
        redispatch_to "/view/";
      },

      ## In apache /viewajax maps to /cachefill
      #  because we want generate the view synchronously to the request
      #  and fill in memcached after its generated
      sub (/cachefill/...) {
        $_[1]->{AJAX_REFRESH} = 2;
        redispatch_psgi($app{'Cache.psgi'}, $_[1]);
      },

      ## This is triggered as an ajax request from the cache-miss page
      sub (/cachetrigger/...) {
        $_[1]->{AJAX_REFRESH} = 1;
        redispatch_psgi($app{'Cache.psgi'}, $_[1]);
      },

      ## In apache /view maps to /cache which will show the cache-miss
      #  page if necessary.
      sub (/cache/...) {
        redispatch_psgi ($app{'Cache.psgi'}, $_[1]);
      },

      sub (/viewajax/...) {
        redispatch_psgi($app{'Rest.psgi'}, $_[1]);
      },

      ($always_memcache ? (
      sub (/view/...) {
        redispatch_psgi($app{'Cache.psgi'}, $_[1]);
      },
      ) : (
      ## this exists so the embedded web server can run without caching
      sub (/view/...) {
        redispatch_psgi ($app{'Rest.psgi'}, $_[1]);
      },
      )),

      ## dump the psgi environment, for testing
      sub (/dump/...) {
        redispatch_psgi ($app{'Dump.psgi'}, $_[1]);
      },

      ## send the browser to the finder view of Genome
      sub (/) {
        redirect_to("/view/genome/search/status.html", $_[1]);
      }
};


Genome::Model::Command::Services::WebApp::Main->run_if_script;
