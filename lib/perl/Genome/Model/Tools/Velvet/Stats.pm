package Genome::Model::Tools::Velvet::Stats;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Velvet::Stats {};

sub help_brief {
    'DEPRECATED - use "gmt velvet metrics"';
}

sub help_detail {
    return;
}

sub execute {
    my $self = shift;
    $self->error_message('Stats is deprecated and will be removed. Please use "gmt velvet metrics" for the same functionality');
    return ;
}

1;

