package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::IndexIllumina; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::IndexIllumina {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes',
};

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

sub _after_update_value {
    my $self = shift;

    my $genome_property_name = $self->genome_property_name;
    if ( $genome_property_name eq 'lane' or $genome_property_name eq 'index_sequence' ) {
        my $update_subset_name = $self->_update_subset_name;
        return if not $update_subset_name;
    }

    return 1;
}

sub _update_subset_name {
    my $self = shift;

    my $genome_entity = $self->genome_entity;
    my $lane = $genome_entity->lane;
    return 1 if not defined $lane;

    my $index_sequence = $genome_entity->index_sequence;
    my $new_subset_name = ( defined $index_sequence )
    ? $lane.'-'.$index_sequence
    : $lane;
    $genome_entity->subset_name($new_subset_name);

    return 1;
}

1;

