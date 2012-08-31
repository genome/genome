package Genome::DruggableGene::View::SearchGo::Html;

use strict;
use warnings;
require UR;

class Genome::DruggableGene::View::SearchGo::Html {
    is => 'Genome::View::Static::Html',
    has_constant => [
        toolkit     => { value => 'html' },
        perspective => { value => 'search-go' },
        desired_perspective => { value => 'status' },
    ],
};

sub _title {
    return 'The Drug-Gene Interactions Database Potentially Druggable Gene Families';
}

sub _html_body {
    return Genome::Sys->read_file(__FILE__ .'.html');
}

1;
