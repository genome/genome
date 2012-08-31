package Genome::Model::Tools::Germline;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::Germline {
    is => 'Command',
    has => [
   ],
};

sub help_brief {
    "Germline Tool",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools germline...    
EOS
}

sub help_detail {
    return <<EOS
Germline!!!
EOS
}


1;

