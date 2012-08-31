package Genome::Command::AddRemoveProperty;

use strict;
use warnings;

use Genome;
      
use Data::Dumper 'Dumper';

class Genome::Command::AddRemoveProperty {
    is => 'Genome::Command::Base',
    is_abstract => 1,
    doc => 'CRUD update property command class.',
};

sub _add_or_remove { Carp::confess('Please use CRUD or implement _add_or_remove in '.$_[0]->class); }
sub _to_or_from { return $_[0]->_add_or_remove eq 'add' ? 'to' : 'from'; }
sub _target_name { Carp::confess('Please use CRUD or implement _target_name in '.$_[0]->class); }
sub _target_name_pl { Carp::confess('Please use CRUD or implement _target_name_pl in '.$_[0]->class); }
sub _target_name_pl_ub { Carp::confess('Please use CRUD or implement _target_name_pl_ub in '.$_[0]->class); }
sub _property_name { Carp::confess('Please use CRUD or implement _property_name in '.$_[0]->class); }
sub _property_name_pl { Carp::confess('Please use CRUD or implement _property_name_pl in '.$_[0]->class); }

sub help_brief {
    return join(' ', split('_', $_[0]->_property_name), $_[0]->_to_or_from, $_[0]->_target_name_pl);
}

sub help_detail {
    my $class = shift;
    return ucfirst($class->_add_or_remove).' '.$class->_property_name.' on '.$class->_target_name_pl.". Parameters are resolved via string.\n";
}

sub execute { # ub = under bar; pl = plural;
    my $self = shift;

    my $target_name = $self->_target_name;
    my $target_name_pl_ub = $self->_target_name_pl_ub;
    my @objects = $self->$target_name_pl_ub;
    my $property_name = $self->_property_name;
    my $property_name_pl = $self->_property_name_pl;
    my $add_or_remove = $self->_add_or_remove;
    my $method = $add_or_remove.'_'.$property_name;
    my $to_or_from = $self->_to_or_from;

    $self->status_message(ucfirst($add_or_remove)." $property_name_pl $to_or_from $target_name_pl_ub...");
    my @new_values = $self->values;
    OBJECT: for my $object ( @objects ) {
        my $object_name = $self->_display_name_for_value($object);
        my @old_values = $object->$property_name_pl;
        my $old_value_names = $self->_display_name_for_value(\@old_values);
        for my $new_value ( @new_values ) {
            my $new_value_name = $self->_display_name_for_value($new_value);
            $self->status_message("$object_name $add_or_remove $property_name '$new_value_name'");
            $object->$method($new_value);
        }
    }
    $self->status_message('Update complete. Committing...');

    return 1; 
}

1;

