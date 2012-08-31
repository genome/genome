package Genome::ModelGroup::View::Metrics::Html;

use strict;
use warnings;
require UR;

class Genome::ModelGroup::View::Metrics::Html {
    is => 'Genome::View::Static::Html',
    has_constant => [
        toolkit     => { value => 'html' },
        perspective => { value => 'metrics' },
        desired_perspective => { value => 'status' },
    ],
};

sub _title {
    return 'Build Metrics For This Model Group';
}

sub _html_body {
    return Genome::Sys->read_file(__FILE__ .'.html');
}

1;

