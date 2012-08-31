package Genome::View::Home::Html;

use strict;
use warnings;
require UR;

class Genome::View::Home::Html {
    is => 'Genome::View::Static::Html',
    has_constant => [
        toolkit     => { value => 'html' },
        perspective => { value => 'home' },
        desired_perspective => { value => 'status' },
    ],
};

sub _title {
    return 'Home';
}

sub _html_body {
    return Genome::Sys->read_file(__FILE__ .'.html');
}

1;
