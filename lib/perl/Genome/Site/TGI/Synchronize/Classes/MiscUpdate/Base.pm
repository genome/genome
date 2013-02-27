package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Base; 

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Base {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate',
};

sub perform_update {
    my $self = shift;

    # Validate
    my $genome_property_name = $self->genome_property_name;
    return $self->failure if not $genome_property_name;

    my $lims_table_name = $self->lims_table_name;
    return $self->failure if not $lims_table_name;

    my $site_tgi_class_name = $self->site_tgi_class_name;
    return $self->failure if not $site_tgi_class_name;

    my $genome_entity = $self->genome_entity;
    return $self->failure if not $genome_entity;

    my $current_value = $self->_get_current_value; # get and set for errors
    $self->{_current_value} = $current_value;

    if ( not grep { $genome_property_name eq $_ } $site_tgi_class_name->properties_to_keep_updated ) {
        $self->status_message('Update for genome property name not supported => '.$genome_property_name);
        return $self->skip;
    }

    # Get values
    my $new_value = $self->new_value;
    my $old_value = $self->old_value;

    # NEW and CURRENT are NULL
    if ( not defined $new_value and not defined $current_value ) {
        $self->is_reconciled(1);
        return $self->success;
    }

    # NEW and CURRENT are NOT NULL and the same
    if ( defined $current_value and defined $new_value and $current_value eq $new_value ) {
        $self->is_reconciled(1);
        return $self->success;
    }

    # OLD and CURRENT value do not match
    if ( defined $old_value and defined $current_value and $old_value ne $current_value ) {
        $self->error_message("Current APipe value ($current_value) does not match the LIMS old value ($old_value)!");
        return $self->failure;
    }

    # Update
    my $updated_value = $self->_update_value;
    return $self->failure if not defined $updated_value;

    if ( not $updated_value or $updated_value ne $new_value ) {
        $self->error_message('Failed to set new value!');
        return $self->failure;
    }

    # YAY!
    return $self->success;
}

sub _get_current_value {
    my $self = shift;

    my $genome_property_name = $self->genome_property_name;
    return if not $genome_property_name;

    my $current_attr = $self->_get_attribute($genome_property_name);
    return $current_attr->attribute_value if $current_attr;

    return;
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

1;

