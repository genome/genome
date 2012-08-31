package Genome::View::SampleManagement::Html;

use strict;
use warnings;
require UR;

class Genome::View::SampleManagement::Html {
    is => 'Genome::View::Static::Html',
    has_constant => [
        toolkit     => { value => 'html' },
        perspective => { value => 'home' },
        desired_perspective => { value => 'status' },
    ],
};

sub _title {
    return 'Sample Management';
}

sub _html_body {
    return Genome::Sys->read_file(__FILE__ .'.html');
}

1;
