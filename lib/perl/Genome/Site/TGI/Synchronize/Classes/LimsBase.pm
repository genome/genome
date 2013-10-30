package Genome::Site::TGI::Synchronize::Classes::LimsBase; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsBase { 
    is => 'UR::Object',
};

sub genome_class_for_comparison {
    my $self = shift;
    return $self->genome_class_for_create;
}

sub genome_class_for_create {
    my $self = shift;
    return 'Genome::'.Genome::Utility::Text::string_to_camel_case($self->entity_name);
}

sub create_in_genome {
    my $self = shift;

    my %params;
    for my $name ( $self->properties_to_copy ) {
        my $value = $self->$name;
        next if not defined $value;
        $params{$name} = $value;
    }

    my $genome_class = $self->genome_class_for_create;
    my $genome_object = eval { $genome_class->create(%params); };
    Carp::confess("Could not create new object of type $genome_class based on object of type " .
    $self->class . " with id " . $self->id . ":\n$@") unless $genome_object;

    return $genome_object;
}

1;

