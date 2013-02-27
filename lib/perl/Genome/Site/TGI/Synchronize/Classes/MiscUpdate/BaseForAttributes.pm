package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Base',
};

sub _get_current_value {
    my $self = shift;

    my $genome_entity = $self->genome_entity;
    return if not $genome_entity;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    my $current_attr = $self->_get_attribute($genome_property_name);
    return if not $current_attr;
    
    return $current_attr->attribute_value;
}

sub _update_value {
    my $self = shift;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    my $current_attr = $self->_get_attribute($genome_property_name);
    $current_attr->delete if $current_attr;

    my $new_attr = $self->_update_attribute($genome_property_name, $self->new_value);
    return if not $new_attr;

    if ( my $after_update_value = $self->can('_after_update_value') ) {
        my $rv = $after_update_value->($self);
        return if not $rv;
    }

    return $new_attr->attribute_value;
}

sub _update_attribute {
    my ($self, $attribute_label, $new_value) = @_;

    my $current_attr = $self->_get_attribute($attribute_label);
    if ( $current_attr ) {
        if ( $current_attr->attribute_value eq $new_value ) {
            return $current_attr;
        }
        $current_attr->delete;
    }

    my $new_attr = $self->_create_attribute($attribute_label, $new_value);
    return $new_attr;
}


sub _get_attribute {
    my ($self, $attribute_label) = @_;

    return Genome::SubjectAttribute->get(
        subject_id => $self->subject_id,
        attribute_label => $attribute_label,
        nomenclature => 'WUGC',
    );
}

sub _create_attribute {
    my ($self, $attribute_label, $attribute_value) = @_;

    my $new_attr = Genome::SubjectAttribute->create(
        subject_id => $self->subject_id,
        attribute_label => $attribute_label,
        attribute_value => $attribute_value,
        nomenclature => 'WUGC',
    );
    if ( not $new_attr ) {
        $self->error_message('Failed to create instrument data attribute!');
        return;
    }

    return $new_attr;
}

1;

