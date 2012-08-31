#!/usr/bin/env genome-perl

use Web::Simple 'HelloWorld';

package HelloWorld;

use Data::Dumper;

dispatch {

    sub (/**) {
        [
            200,
            [ 'Content-type', 'text/plain' ],
            [ Data::Dumper->new( [ \@_ ] )->Dump ]
        ];
      }, sub (GET) {
        [ 200, [ 'Content-type', 'text/plain' ], ['hello hello!'] ];
      }, sub () {
        [ 405, [ 'Content-type', 'text/plain' ], ['Method not allowed'] ];
      }
};

HelloWorld->run_if_script;
