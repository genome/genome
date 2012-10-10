#!/usr/bin/env genome-perl

use Web::Simple 'Genome::Model::Command::Services::WebApp::Rest';

package Genome::Model::Command::Services::WebApp::Rest;

#my $res_path = Genome::Model::Command::Services::WebApp->res_path;

our $loaded = 0;
sub load_modules {
    return if $loaded;
    eval "
        use above 'Genome';
        use Plack::MIME;
        use Plack::Util;
        use Plack::Request;
        use Cwd;
        use HTTP::Date;
        use JSON;
        use UR::Object::View::Default::Xsl qw/type_to_url url_to_type/;
    ";
    if ($@) {
        die "failed to load required modules: $@";
    }

    # search's callbacks are expensive, web server can't change anything anyway so don't waste the time
    Genome::Search->unregister_callbacks();
}



my $new_print = sub { 

# this attempt to make it do a smarter "print" doesn't work so just make a no-op right now
=cut
    my $self = shift; 
    $DB::single = 1;
    # 6 = FCGI::Stream::FCGI_STDOUT constant
    if ($self->{type} eq 6) {
        STDERR->print("BAD STDOUT:", @_);
    }
    else {
        $old_print->($self,@_);
    }
=cut
};

sub dispatch_request {

    sub (DELETE + /** + %@*) {

        load_modules();
        my ($self, $url, $params) = @_;
        my $class = url_to_type($url);
        my @ids = $params->{'ids[]'};

        my $status;
        eval {
            local *FCGI::Stream::PRINT = $new_print;
            my @objs;
            if (@ids > 1 || (@ids && $ids[0] != undef)) {
                @objs = $class->get(@ids);
            } else {
                my $obj = $class->get($params);
                @objs = ($obj);
            }
            if (@objs > 200) { die 'remove this when youre ready to delete lots of stuff from a url ;)'; }
            for my $obj (@objs) {
                $obj->delete();
            }
            $status = UR::Context->commit();
        };

        my ($code, $body);
        my $r = {};
        if ($@ || !$status) {
            $r->{'error'} = substr(UR::Context->error_message() || $@,0,5000);
            $code = 200;
            $body = to_json( $params, { ascii => 1, allow_nonref => 1, });
        } else { 
            # worked- later this could be 202 for asyncronous (via task) deletion
            $code = 200;
            $r->{'deleted_class'} = $class;
            $r->{'deleted_ids'} = \@ids;
            $body = to_json( $r, { ascii => 1, allow_nonref => 1, });
        }

        return [$code, ['Content-type' => 'text/plain'], [$body]];
    }, 
    sub (PUT + /** + %*) {

        load_modules();

        my ($self, $url, $params) = @_;
        my ($code, $obj);
        my $class = url_to_type($url); # UR::Object::View::Default::Xsl

        my $status;
        eval { 
            local *FCGI::Stream::PRINT = $new_print;
            $obj = $class->create(%$params);
            $status = UR::Context->commit();
        };
        my $error = $@;

        my $body;
        if ($error || !$status) {
            $params->{'error'} = substr(UR::Context->error_message() || $error,0,5000);
            $code = 200; # OK (didnt work)
            $params->{'code'} = $code;
            $body = to_json( $params, { 
                                ascii => 1,
                                allow_nonref => 1,
                    });
        } else {
            $code = 201; # CREATED
            $params->{'code'} = $code;
            $params->{'id'} = $obj->id();
            $body = to_json( $params, { 
                                ascii => 1,
                                allow_nonref => 1,
                    });

        }

#        my $size; { use bytes; $size = length($body); }
#                'Content-Length' => $size

        return [
            $code,
            [
                'Content-type'   => "text/plain"
            ],
            [$body]
        ];

    },

    sub (POST + /** + %*) {
        load_modules();
        my ($self, $url, $params) = @_;
        my ($code, $obj, $body);
        my $class = url_to_type($url); # UR::Object::View::Default::Xsl

        my $id = $params->{'id'};
        $obj = $class->get($id);
        if (!$id) {
            $code = 400;
            $body = "ERROR: No ID passed in.";
            return [ $code, [ 'Content-type'   => "text/plain" ], [$body] ];
        }
    
        if (!$obj) {
            $code = 404;
            $body = "ERROR: No object by that ID found.";
            return [ $code, [ 'Content-type'   => "text/plain" ], [$body] ];
        }
        eval {
            for my $p (keys %$params) {
                 $obj->$p($params->{$p});
             }
            UR::Context->commit;
        };

        if ($@) {
            $code = 500; 
            return [ $code, [ 'Content-type'   => "text/plain" ], ["An error occurred processing the input: $@"] ];
        } else {
            $code = 201; # CREATED
            $params->{'id'} = $obj->id();
        }

        $body = to_json( $params, {
                                ascii => 1,
                                allow_nonref => 1,
                    });
        return [
            $code,
            [
                'Content-type'   => "text/plain"
            ],
            [$body]
        ];

    },

    # Matcher for Static content related to a view
    # **/ = class name
    # */ = perspective.toolkit
    # * + .* = filename & extension
    # matches urls like /view/Genome/Model/status.html/foo.jpg
    # would map to Genome/Model/View/Status/Html/foo.jpg

    sub (GET + /**/*/* + .*) {
        # these get passed in from the matcher as documented above!
        my ( $self, $class, $perspective_toolkit, $filename, $extension ) = @_;

        $filename = $filename . "." . $extension;

        load_modules();

        if ( $class =~ /\./ ) {

            # matched on some multi-part path after the perspective_toolkit
            if ( $class =~ s/\/(.+?)$//g ) {
                $filename            = $perspective_toolkit . '/' . $filename;
                $perspective_toolkit = $1;

                if ( $perspective_toolkit =~ s/\/(.+)$//g ) {
                    $filename = $1 . '/' . $filename;
                }
            } else {

                # let some other handler deal with this, i cant parse it
                return;
            }
        } elsif ( index( $perspective_toolkit, '.' ) < 0 ) {

            # doesn't have a period in it?  probably didnt want us to match
            return;
        }

        $class = url_to_type($class);
        my ( $perspective, $toolkit ) = split( /\./, $perspective_toolkit );
        my $mime_type = Plack::MIME->mime_type(".$extension");

        my $view_class = UR::Object::View->_resolve_view_class_for_params(
            subject_class_name => $class,
            perspective        => $perspective,
            toolkit            => $toolkit
        );

        # 404 handler will rewrite text/plain into a prettier format 
        unless ($view_class) {
            return [
                404,
                [ 'Content-type', 'text/plain' ],
                ["No view for $class $perspective $toolkit"]
            ];
        }

        my $base_dir = $view_class->base_dir;
        unless ( -d $base_dir ) {
            return [
                404,
                [ 'Content-type', 'text/plain' ],
                ["No resource directory for $view_class"]
            ];
        }

        my $full_path = $base_dir . '/' . $filename;

        unless ( -e $full_path ) {
            return [
                404,
                [ 'Content-type', 'text/plain' ],
                ["No resource $filename for $view_class"]
            ];
        }

        open my $fh, "<:raw", $full_path
          or return [ 403, [ 'Content-type', 'text/plain' ], ['forbidden'] ];

        my @stat = stat $full_path;

        Plack::Util::set_io_path( $fh, Cwd::realpath($full_path) );

        ## Plack should have set binmode for us, workaround here because it's dumb.
        if ( $ENV{'GATEWAY_INTERFACE'} ) {
            binmode STDOUT;
        }

        return [
            200,
            [
                'Content-type'   => $mime_type,
                'Content-Length' => $stat[7],
                'Last-Modified'  => HTTP::Date::time2str( $stat[9] )
            ],
            $fh
        ];
      },

      # Second matcher maps UR views
      #/** class
      #/* perspective
      #.* toolkit
      # + ?@*   slurp query portion into a hash of arrays (so you could say ?foo=a&foo=b&foo=c and get foo=>[a,b,c])
      sub (GET + /**/* + .* + ?@*) {
        my ( $self, $class, $perspective, $toolkit, $args ) = @_;

        $DB::single = 1;

        load_modules();

        $class = url_to_type($class);
        $perspective =~ s/\.$toolkit$//g;

        my $mime_type = Plack::MIME->mime_type(".$toolkit");

        # flatten these where only one arg came in (don't want x=>['y'], just x=>'y')
        for my $key ( keys %$args ) {
            if ( index( $key, '_' ) == 0 ) {
                delete $args->{$key};
                next;
            }
            my $value = $args->{$key};

            if ( $value and scalar @$value eq 1 ) {
                $args->{$key} = $value->[0];
            }
        }

        my %view_special_args;
        for my $view_key (grep {$_=~ m/^-/} keys %$args) {
            $view_special_args{substr($view_key,1,length($view_key))} = delete $args->{$view_key}; 
        }

        my @matches;

        if ($class->isa("UR::Object::Set")) {
            $class =~ s/::Set$//;
            eval {
                @matches = $class->define_set(%$args);
            };
            if ($@) {
                my $arg_str = join(' ', %$args);
                return [ 500, [ 'Content-type', 'text/plain' ],
                    ["Failed to define_set() for $class with args: $arg_str"] ];
            }
        }
        else {
            eval {
                @matches = $class->get(%$args);
            };
            if ($@) {
                my $arg_str = join(' ', %$args);
                return [ 500, [ 'Content-type', 'text/plain' ],
                    ["Failed to get() for $class with args: $arg_str"] ];
            }
        }

        unless (@matches) {
            return [ 404, [ 'Content-type', 'text/plain' ],
                ['No object found'] ];
        }
        die 'matched too many; list not yet supported' unless ( @matches == 1 );

        my %view_args = (
            perspective => $perspective,
            toolkit     => $toolkit
        );

        if ( $toolkit eq 'xsl' || $toolkit eq 'html' ) {
            $view_args{'xsl_root'} =
              Genome->base_dir . '/xsl';    ## maybe move this to $res_path?
            $view_args{'xsl_path'} = '/static/xsl';
            $view_args{'html_root'} =
              Genome->base_dir . '/View/Resource/Html/html';


            #            $view_args{'rest_variable'} = '/view';

            $view_args{'xsl_variables'} = {
                rest      => '/view',
                resources => '/view/genome/resource.html'
            };
        }

        my $view;
        # all objects in UR have create_view
        # this probably ought to be revisited for performance reasons because it has to do a lot of hierarchy walking
        eval { $view = $matches[0]->create_view(%view_args, %view_special_args); };

        if ( $@ && !$view ) {
            $view_args{'desired_perspective'} = $perspective;
            $view_args{'perspective'}         = 'default';

            eval { $view = $matches[0]->create_view(%view_args); };
            if ($@) {
                return [
                    404, [ 'Content-type', 'text/plain' ],
                    ['No view found']
                ];
            }
        }

        die 'no_view' unless ($view);

        my $content = $view->content();
        
        # the 'cacheon' cookie is handled by the javascript /View/Resource/Html/js/app/ui-init.js to tell the browser it's a
        # cached page and show the timer.  
        #this cookie has an expiration date in the past, so it tells the browser to expire it right now...
        [ 200, [ 'Content-type', $mime_type, 'Set-Cookie', 'cacheon=1; expires=Thu, 01-Jan-2000 00:00:00 GMT' ], [$content] ];
      }
};

Genome::Model::Command::Services::WebApp::Rest->run_if_script;
