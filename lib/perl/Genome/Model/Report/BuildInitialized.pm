package Genome::Model::Report::BuildInitialized;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Report::BuildInitialized {
    is => 'Genome::Model::Report::BuildEventBase',
};

sub _add_to_report_xml {
    my $self = shift;

    $self->_add_build_event_dataset
        or return;

    return 1;
}

1;

#$HeadURL$
#$Id$
