package Genome::Model::Tools::Dgidb::Base;

use strict;
use warnings;
use Genome;
use LWP::UserAgent;
use HTTP::Request::Common;



my $DOMAIN   = 'http://dgidb.genome.wustl.edu/';

class Genome::Model::Tools::Dgidb::Base {
    is => 'Command',
    is_abstract => 1,
};

sub get_domain {
    return $DOMAIN;
}

sub post_request {
    my $self = shift;
    my $params_ref = shift;
    my $ua = LWP::UserAgent->new;
    return $ua->request(POST $self->get_domain . $self->get_api_path, $params_ref);
}

sub get_api_path {
    my $self = shift;
    $self->error_message("Must override get_api_path in child class");
    return;
}

1;

