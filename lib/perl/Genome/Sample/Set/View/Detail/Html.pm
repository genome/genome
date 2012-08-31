package Genome::Sample::Set::View::Detail::Html;

use strict;
use warnings;
require UR;

class Genome::Sample::Set::View::Detail::Html {
    is => 'Genome::View::Static::Html',
    has_constant => [
        toolkit     => { value => 'html' },
        perspective => { value => 'detail' },
        desired_perspective => { value => 'status' },
    ],
};

sub _title {
    return 'Samples and their Attributes';
}

sub _html_body {
    return Genome::Sys->read_file(__FILE__ .'.html');
}

1;

