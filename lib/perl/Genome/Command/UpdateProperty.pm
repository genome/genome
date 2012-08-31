package Genome::Command::UpdateProperty;

use strict;
use warnings;

use Genome;
      
use Data::Dumper 'Dumper';

class Genome::Command::UpdateProperty {
    is => 'Genome::Command::Base',
    is_abstract => 1,
    doc => 'CRUD update property command class.',
};

sub _target_name_pl { Carp::confess('Please use CRUD or implement _target_name in '.$_[0]->class); }
sub _target_name_pl_ub { Carp::confess('Please use CRUD or implement _target_name_ub in '.$_[0]->class); }
sub _property_name { Carp::confess('Please use CRUD or implement _property_name in '.$_[0]->class); }
sub _property_doc { return "modify '".join(' ', split('_', $_[0]->_property_name))."'"; }
sub _only_if_null { Carp::confess('Please use CRUD or implement _only_if_null in '.$_[0]->class); }

sub help_brief { return $_[0]->_property_doc; }

sub help_detail {
    my $class = shift;
    my $help = 'Update '.$class->_property_name.' on '.$class->_target_name_pl.' Parameters are resolved via string.';
    if ( $class->_only_if_null ) {
        $help .= " This property can only be updated if it is NULL.";
    }
    $help .= "\n\n";
    return $help;
}

sub execute {
    my $self = shift;

    my $property_name = $self->_property_name;
    my $target_name_pl_ub = $self->_target_name_pl_ub;
    my @objects = $self->$target_name_pl_ub;

    my $new_value = $self->value;
    if ( $new_value eq '' ) {
        #FIXME - this will set strings and direct objects to null.
        # It does not work for a property that is through an many property. We are unsure of how we want to handle unsetting values.
        $self->error_message("Cannot set $property_name to NULL. Sorry.");
        return;
    }
    $new_value = undef if $new_value eq '';
    my $new_value_name = $self->_display_name_for_value($new_value);

    $self->status_message("Update $property_name for $target_name_pl_ub...");
    for my $object ( @objects ) {
        my $object_name = $self->_display_name_for_value($object);
        my $old_value = $object->$property_name;
        my $old_value_name = $self->_display_name_for_value($old_value);
        if ( $self->_only_if_null and defined $old_value ) {
            $self->status_message("$object_name  - cannot update because $property_name has value ($old_value_name)!");
            next;
        }
        $object->$property_name($new_value);
        $self->status_message("$object_name - update from '$old_value_name' to '$new_value_name'.");
    }
    $self->status_message('Update complete. Committing changes...');

    return 1; 
}

1;

