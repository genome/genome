package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::RunRegion454; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::RunRegion454 {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Base',
};

class Genome::InstrumentData::RunRegion454 { };
sub Genome::InstrumentData::RunRegion454::get {
    my ($class, $id) = @_;
    my $self = $class->SUPER::get(id => $id);
    return $self if $self;
    return $class->create(id => $id);
}
sub Genome::InstrumentData::RunRegion454::instrument_data_454 {
    my $self = shift;

    my @instrument_data_454 = Genome::InstrumentData::454->get(
        'attributes.attribute_label' => 'region_id',
        'attributes.attribute_value' => $self->id,
        'attributes.nomenclature' => 'WUGC',
    );
    if ( not @instrument_data_454 ) {
        $self->error_message('Failed to find instrument data 454 for region id! '.$self->id);
        return;
    }

    return @instrument_data_454;
};

sub _get_current_value {
    my $self = shift;

    return $self->old_value;
}

sub _update_value {
    my $self = shift;

    my $genome_entity = $self->genome_entity;
    return if not $genome_entity;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    my @instrument_data_454 = $genome_entity->instrument_data_454;
    return if not @instrument_data_454;

    for my $instrument_data_454 ( @instrument_data_454 ) {
        $instrument_data_454->$genome_property_name($self->new_value);
        next if $genome_property_name ne 'region_number';
        $instrument_data_454->$genome_property_name($self->new_value);
        my $region_number = $instrument_data_454->region_number;
        next if not defined $region_number;
        my $index_sequence = $instrument_data_454->index_sequence;
        my $new_subset_name = ( defined $index_sequence )
        ? $region_number.'-'.$index_sequence
        : $region_number;
        $instrument_data_454->subset_name($new_subset_name);
    }

    if ( my $after_update_value = $self->can('_after_update_value') ) {
        my $rv = $after_update_value->($self);
        return if not $rv;
    }

    return 1;
}

sub _validate_value_set {
    return 1;
}

1;

