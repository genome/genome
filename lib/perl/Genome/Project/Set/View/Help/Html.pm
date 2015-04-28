package Genome::Project::Set::View::Help::Html;

use strict;
use warnings;
require UR;

class Genome::Project::Set::View::Help::Html {
    is => 'Genome::View::Static::Html',
    has_constant => [
        toolkit     => { value => 'html' },
        perspective => { value => 'help' },
        desired_perspective => { value => 'status' },
    ],
};

sub _title {
    return 'How to use projects';
}

sub _html_body {
    my $body = Genome::Sys->read_file(__FILE__ .'.html');

    # Crappy templating
    my $email = Genome::Config::get('email_pipeline');
    $body =~ s/##EMAIL##/$email/g;
    return $body;
}

1;

