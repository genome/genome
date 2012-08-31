package Genome::Model::Event::Build::RnaSeq;

#REVIEW fdu 11/19/2009
#command_subclassing_model_property becomes obsolete because its
#subclass implement with their own command_subclassing_model_property

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::RnaSeq {
    is => 'Genome::Model::Event',
};

sub command_subclassing_model_property {
    return 'sequencing_platform';
}

1;

