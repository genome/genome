package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Base',
};

sub _get_current_attr {
    my $self = shift;

    my $genome_entity = $self->genome_entity;
    return if not $genome_entity;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    return $genome_entity->attributes(
        attribute_label => $genome_property_name,
        nomenclature => 'WUGC',
    );
}

sub _get_current_value {
    my $self = shift;

    my $genome_entity = $self->genome_entity;
    return if not $genome_entity;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    my $current_value = eval{ $genome_entity->$genome_property_name; };
    return $current_value if defined $current_value;

    my $current_attr = $self->_get_current_attr;
    return if not $current_attr;
    
    return $current_attr->attribute_value;
}

sub _update_value {
    my $self = shift;

    my $genome_entity = $self->genome_entity;
    return if not $genome_entity;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;
    my $pmeta = $genome_entity->__meta__->property($genome_property_name);
    unless (defined $pmeta->{via} and $pmeta->{via} eq "attributes") {
        return $self->SUPER::_update_value();
    }

    my $current_attr = $self->_get_current_attr;
    $current_attr->delete if $current_attr;

    my $new_attr = $genome_entity->add_attribute(
        attribute_label => $genome_property_name,
        attribute_value => $self->new_value,
        nomenclature => 'WUGC',
    );
    return if not $new_attr;

    if ( my $after_update_value = $self->can('_after_update_value') ) {
        my $rv = $after_update_value->($self);
        return if not $rv;
    }

    return $new_attr->attribute_value;
}

1;

