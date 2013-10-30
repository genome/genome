package Genome::Site::TGI::Synchronize::Classes::LimsInstDataBase; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsInstDataBase { 
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
};

sub genome_class_for_create {
    my $self = shift;
    return 'Genome::InstrumentData::'.Genome::Utility::Text::string_to_camel_case($self->entity_name);
}

1;

