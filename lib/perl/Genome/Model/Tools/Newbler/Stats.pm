package Genome::Model::Tools::Newbler::Stats;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Newbler::Stats {};

sub help_brief {
    'DEPRECATED - use "gmt newbler metrics"';
}

sub help_detail {
    return;
}

sub execute {
    my $self = shift;
    $self->error_message('Stats is deprecated and will be removed. Please use "gmt newbler metrics" for the same functionality');
    return ;
}

1;

