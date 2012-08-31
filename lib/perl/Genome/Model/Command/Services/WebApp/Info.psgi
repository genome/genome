#!/usr/bin/env genome-perl

use Web::Simple 'Genome::Model::Command::Services::WebApp::Info';


package Genome::Model::Command::Services::WebApp::Info;

use English;
use Data::Dumper;

sub dispatch_request {
    sub (/**) {
        [
            200,
            [ 'Content-type', 'text/plain' ],
            [ dump_info() ]
        ];
      }, sub () {
        [ 405, [ 'Content-type', 'text/plain' ], ['Method not allowed'] ];
      }
};

Genome::Model::Command::Services::WebApp::Info->run_if_script;


sub dump_info {

    my @inc = map {join(': ',$_,$INC{$_})} sort keys %INC;
    my @env = map {join(': ',$_,$ENV{$_})} sort keys %ENV;

    my $i = join("\n",
        'This page dumps %INC and %ENV for debugging.',
        "\n",
        scalar(localtime()),
        "EXECUTABLE_NAME: " . $^X,
        "PROGRAM NAME: $PROGRAM_NAME",
        "\n***INC***",
        join("\n", @inc),
        "\n\n***ENV***",
        join("\n", @env)
    );

    return $i;
}


