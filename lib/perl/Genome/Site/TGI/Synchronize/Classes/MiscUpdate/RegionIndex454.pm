package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::RegionIndex454; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::RegionIndex454 {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Base',
};

sub _get_current_value {
    my $self = shift;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    my $current_attr = $self->_get_attribute($genome_property_name);
    return $current_attr->attribute_value if $current_attr;

    return;
}

sub _get_attribute {
    my ($self, $attribute_label) = @_;

    return Genome::InstrumentDataAttribute->get(
        instrument_data_id => $self->subject_id,
        attribute_label => $attribute_label,
        nomenclature => 'WUGC',
    );
}

sub _create_attribute {
    my ($self, $attribute_label, $attribute_value) = @_;

    my $new_attr = Genome::InstrumentDataAttribute->create(
        instrument_data_id => $self->subject_id,
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

sub _update_value {
    my $self = shift;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    my $current_attr = $self->_get_attribute($genome_property_name);
    $current_attr->delete if $current_attr;

    my $new_attr = $self->_update_attribute($genome_property_name, $self->new_value);
    return if not $new_attr;

    if ( $genome_property_name eq 'region_number' or $genome_property_name eq 'index_sequence' ) {
        my $update_subset_name = $self->_update_subset_name;
        return if not $update_subset_name;
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

sub _update_subset_name {
    my $self = shift;
    $DB::single=1;

    my $genome_entity = $self->genome_entity;
    my $region_number = $genome_entity->region_number;
    return 1 if not defined $region_number; # should not happen

    my $index_sequence = $genome_entity->index_sequence;
    my $new_subset_name = ( defined $index_sequence )
    ? $region_number.'-'.$index_sequence
    : $region_number;
    $genome_entity->subset_name($new_subset_name);

    return 1;
}

1;

