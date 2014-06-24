#!/usr/bin/env genome-perl
package Genome::Model::Command::Services::WebApp::Cache;

use strict;
use warnings;
use Data::Dumper;

# don't cache static things.
# every url gets matched against this, so it shouldn't be a long list
# if it becomes long, consider a redesign.

#our @never_cache = (
#    qr{(?<!html)$},
#    qr{Genome::Search::Query}i,
#    qr{genome/search}i
#);

our @always_cache = (
    qr{genome/model}i,
    qr{genome/sample/}i
);

use Plack::Util;
use Digest::MD5;
use Cache::Memcached;
use Storable qw/freeze thaw/;
use Sys::Hostname qw/hostname/;

our $environment = (hostname eq 'vm44' || hostname eq 'vm62.gsc.wustl.edu') ? 'prod' : 'dev';
our %servers = ('prod' => $ENV{'GENOME_SYS_SERVICES_MEMCACHE'}, 
                'dev' => $ENV{'GENOME_SYS_SERVICES_MEMCACHE'}, 
                'local' => 'localhost:11211'
                );
our $cache_timeout = 0;
our $lock_timeout = 600;
our $server;

sub environment {
    my $class = shift;
    $environment = shift;
    undef $server;
}

sub server {
    my $class = shift;

    unless (defined $server) {
        $server = Cache::Memcached->new({
            servers => [$servers{$environment}],
            debug => 0,
            compress_threshold => 10_000
        });
    }
    return $server;
}

sub hash_url {
    my $class = shift;
    my $url = shift;

    return Digest::MD5::md5_base64($url);
}

sub cache_key_for_url {
    my $class = shift;
    my $url = shift;

    return 'genome_wac:' . $class->hash_url($url);
}

sub lock_key_for_url {
    my $class = shift;
    my $url = shift;

    return 'genome_lock:' . $class->hash_url($url);
}

sub get {
    my $class = shift;
    my $url = shift;

    return $class->server->get($class->cache_key_for_url($url));
}

sub set {
    my $class = shift;
    my $url = shift;
    my $value = shift;

    return $class->server->set($class->cache_key_for_url($url),$value,$cache_timeout);
}

sub lock {
    my $class = shift;
    my $url = shift;

    return $class->server->add($class->lock_key_for_url($url),$$,$lock_timeout);
}

sub unlock {
    my $class = shift;
    my $url = shift;

    return $class->server->delete($class->lock_key_for_url($url));
}

sub getlock {
    my $class = shift;
    my $url = shift;

    return $class->server->get($class->lock_key_for_url($url));
}

sub delete {
    my $class = shift;
    my $url = shift;

    return $class->server->delete($class->cache_key_for_url($url));
}


sub {
    my $class = __PACKAGE__;
    my ($env) = @_;

    my $ajax_refresh = $env->{AJAX_REFRESH};

    my $url = $env->{'PATH_INFO'};
    if (exists $env->{'QUERY_STRING'} && defined $env->{'QUERY_STRING'}) {
        $url .= '?' . $env->{'QUERY_STRING'};
    }

    my $gen = sub {
        my $rest_app = $Genome::Model::Command::Services::WebApp::Main::app{'Rest.psgi'};
        my $resp;
        if ($class->lock($url)) {

            # Looks like this URL isnt in cache yet- cache it

            ## override HTTP_ACCEPT to tell it we want html
            $env->{HTTP_ACCEPT} = "application/xml,application/xhtml+xml,text/html";

            $resp = Plack::Util::run_app $rest_app, $env;

            if ($resp->[0] == 500) {
                $class->unlock($url);
                return $resp;
            }
            if ( ref($resp->[2]) eq 'ARRAY') {
                my $found = 0;
                for (my $i=0; $i < scalar(@{ $resp->[1] }); $i += 2) {
                    if ($resp->[1][$i] eq 'Set-Cookie') {
                        $resp->[1][$i+1] = 'cacheon=1';
                        $found=1;
                        last;
                    }
                }
                if (!$found) {
                    push @{ $resp->[1] }, 'Set-Cookie' => 'cacheon=1';
                }

                $resp->[3] = time;
                if (!$class->set($url,freeze($resp))) {
                    $class->unlock($url);

                    return [
                        500,
                        [ 'Content-type' => 'text/html' ],
                        [ 'Memcached is down' ]
                    ];
                }
            }
            $class->unlock($url);
        } else {

            # Looks like this URL is already in cache

            my $v;
            do {
                $v = $class->getlock($url);
                sleep 1 if $v;
            } while ($v);

            if (defined wantarray) {
                my $v = $class->get($url);
                $resp = thaw($v);
            }
        }

        return $resp;
    };

    if (defined $ajax_refresh && $ajax_refresh == 1) {

        my $resp = $gen->();
        if ($resp->[0] == 500) {
            return $resp; 
        }

        # don't send back the content because we don't care, just want to tell the caller it's ready to ask for
        return [
            200,
            [ 'Content-type' => 'text/html' ],
            [ 'Done' ]
        ];
    } elsif (defined $ajax_refresh && $ajax_refresh == 2) {
        ## ajax request wants a to wait for the page to be generated
        #  without a placeholder to placate the user

        my $resp;
        if (my $v = $class->get($url)) {

            my $no_cache = $env->{'HTTP_CACHE_CONTROL'} || $env->{'HTTP_PRAGMA'};
            # opera always sends no-cache header, so we can't reliably detect
            # a shift-reload
            if ($env->{'HTTP_USER_AGENT'} =~ /Opera/) {
                $no_cache = '';
            }
            
            # this stuff is due to Opera's freakiness as well
            if (exists $env->{'HTTP_X_MAX_AGE'}) {
                my $age = $env->{'HTTP_X_MAX_AGE'};

                $resp = thaw($v);

                if (!$resp->[3] || (time - $resp->[3]) > $age) {
                    undef $resp;
                }

            } elsif (defined $v && $no_cache ne 'no-cache') {
                $resp = thaw($v);
            }
        }

        if (!$resp) {
            $resp = $gen->();
        }

        return [@$resp[0,1,2]];
    } else {
        my $skip_cache = 1;
        
        for my $re (@always_cache) {
            if ($env->{'PATH_INFO'} =~ $re) {
                $skip_cache = 0;
                last;
            }
        }

        if ($skip_cache) {
            my $rest_app = $Genome::Model::Command::Services::WebApp::Main::app{'Rest.psgi'};
            my $resp = Plack::Util::run_app $rest_app, $env;

            return $resp;
        }

        my $v = $class->get($url);

        my $no_cache = $env->{'HTTP_CACHE_CONTROL'} || $env->{'HTTP_PRAGMA'};

        # opera always sends no-cache header, so we can't reliably detect
        # a shift-reload
        if ($env->{'HTTP_USER_AGENT'} =~ /Opera/) {
            $no_cache = '';
        }

        # if we don't have a no-cache and it got pulled from memcache, thaw it and hand it back
        # otherwise spit back a page with an ajax call back to test if we have content that's ready
        # (i.e. call to /cachetrigger)
        my $s = eval { thaw($v) }; # eval catches thaw error due to Perl 5.8 to 5.10 migration
        if (defined $v && defined $no_cache && $no_cache ne 'no-cache' && defined $s) {
            return [@$s[0,1,2]];
        } else {
            my $content = q[
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
    <!--template: status/root.xsl:match "/"-->
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
    <title>Caching Page</title>
    <link rel="shortcut icon" href="/res/img/gc_favicon.png" type="image/png" />
    <link rel="stylesheet" href="/res/css/blueprint/screen.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/blueprint/print.css" type="text/css" media="print" />
    <link rel="stylesheet" href="/res/css/master.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/buttons.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/icons.css" type="text/css" media="screen, projection" />
    <link rel="stylesheet" href="/res/css/forms.css" type="text/css" media="screen, projection" />
    <link type="text/css" href="/res/css/jquery-ui.css" rel="stylesheet" />
    <link href="/res/css/jquery-ui-overrides.css" type="text/css" rel="stylesheet" media="screen, projection" />
  <script type="text/javascript" src="/res/js/pkg/jquery.js"></script>

  <script type="text/javascript">
   (function($) {

     $(document).ready(function() {

        $("#ajax_status")
        .addClass('success')
        .bind("ajaxSend", function(){
            $(this).removeClass('success error').addClass('loading').html('Loading').show();
        })
        .bind("ajaxSuccess", function(){
            $(this).removeClass('loading').addClass('success').html('Success').hide('slow');
        })
        .bind("ajaxError", function(){
            var mailtoAddr = 'apipe-support@rt.gsc.wustl.edu';
            var mailtoSubj = 'Error Loading Page';
            var mailtoBody = 'URL: ' + encodeURIComponent(window.location.href);
            var mailtoLink = encodeURI('mailto:' + mailtoAddr + '?subject=' + mailtoSubj + '&body=' + mailtoBody);
            var theError   = 'Error loading page. Please email <a href="' + mailtoLink + '">' + mailtoAddr + '</a> if you believe this page should exist.';
            $(this).removeClass('loading').addClass('error').html(theError);
        })
        .hide();

       $.ajax({
         url: '/cachetrigger] . $url . q[',
         success: function(data) {
           location.reload();
         }
       });
     });

   })(jQuery)
  </script>
 </head>
 <body>
  <div class="page" style="width: 500px;padding-top: 45px;">
    <div class="content rounded shadow" style="padding-top: 0;" >
      <div class="header rounded-top gradient-grey">
        <div class="container" style="width: 480px;">
          <div class="title app_cache_miss_32">
            <h1>Loading View</h1>
          </div>
        </div>
      </div>

    <div class="container" style="width: 480px;">
      <div class="span-12 last">
        <div class="rounded" style="margin-bottom: 10px;">
          <div class="padding10">
            <p></p>
            <div id="ajax_status"></div>
          </div>
        </div>
      </div>
    </div>
    </div>
  </div>
 </body>
</html>
];

            return [
                200,
                [ 'Content-type' => 'text/html' ],
                [ $content ]
            ];
        }
    }
};
