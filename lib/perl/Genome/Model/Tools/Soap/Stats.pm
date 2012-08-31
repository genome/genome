package Genome::Model::Tools::Soap::Stats;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Soap::Stats {};

sub help_brief {
    'DEPRECATED - use "gmt soap metrics"';
}

sub help_detail {
    return;
}

sub execute {
    my $self = shift;
    $self->error_message('Stats is deprecated and will be removed. Please use "gmt soap metrics" for the same functionality');
    return ;
}

1;

