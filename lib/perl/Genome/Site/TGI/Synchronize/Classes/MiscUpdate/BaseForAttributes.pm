package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::BaseForAttributes {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Base',
};

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

