package Genome::Site::TGI::Synchronize::Classes::LimsInstDataBase; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsInstDataBase { 
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsBase',
};

sub genome_class_for_create {
    my $self = shift;
    my $entity_name = $self->entity_name;
    $entity_name =~ s/instrument data //;
    return 'Genome::InstrumentData::'.Genome::Utility::Text::string_to_camel_case($entity_name);
}

1;

