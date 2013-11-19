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

    my %params = $self->params_for_create_in_genome;
    return if not %params;

    my $genome_class = $self->genome_class_for_create;
    my $genome_object = eval { $genome_class->create(%params); };
    if ( not $genome_object ) {
        Carp::confess("$@\nFailed to create $genome_class with parmas: ".Data::Dumper::Dumper(\%params));
    }

    return $genome_object;
}

sub params_for_create_in_genome {
    my $self = shift;

    my $meta = $self->__meta__;

    my %params;
    for my $name ( $self->properties_to_copy ) {
        my @value = $self->$name;
        next if not @value;
        my $property = $meta->property_meta_for_name($name);
        if ( $property->is_many ) {
            $params{$name} = \@value;
        }
        else {
            $params{$name} = $value[0];
        }
    }

    return %params;
}

sub lims_property_name_to_genome_property_name {
    return $_[1];
}

1;

