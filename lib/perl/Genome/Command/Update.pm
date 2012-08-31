package Genome::Command::Update;

use strict;
use warnings;

use Genome;
      
require Carp;
use Data::Dumper 'Dumper';
require Scalar::Util;

class Genome::Command::Update {
    is => 'Command::V2',
    is_abstract => 1,
    attributes_have => [ is_add_remove => { is => 'Boolean', is_optional => 1, } ],
    doc => 'CRUD update command class.',
};

sub _target_name { Carp::confess('Please use CRUD or implement _target_name in '.$_[0]->class); }
sub _target_name_pl { Carp::confess('Please use CRUD or implement _target_name_pl in '.$_[0]->class); }
sub _target_name_pl_ub { my $target_name_pl = $_[0]->_target_name_pl; $target_name_pl =~ s/ /\_/g; return $target_name_pl }
sub _only_if_null { Carp::confess('Please use CRUD or implement _only_if_null in '.$_[0]->class); }

sub sub_command_sort_position { .3 };

sub help_brief {
    return 'update '.$_[0]->_target_name_pl;
}

sub help_detail {
    my $class = shift;
    my $target_name_pl = $class->_target_name_pl;
    my $help = "This command updates $target_name_pl resolved via text string. Many $target_name_pl can be indicated at one time, as well as many properties to update.\n\n";
    if ( @{$class->_only_if_null} ) {
        $help .= 'These properties can only be updated if NULL: '.join(', ', @{$class->_only_if_null}).'.';
    }
    return $help;
}

sub execute {
    my $self = shift;

    $self->status_message('Update: '.$self->_target_name_pl);

    my $class = $self->class;
    my $target_name_pl_ub = $self->_target_name_pl_ub;
    my $only_if_null = $self->_only_if_null;
    my @objects = $self->$target_name_pl_ub;
    my $properties_requested_to_update = 0;
    my $success = 0;

    PROPERTY: for my $property_meta ( $self->__meta__->property_metas ) {
        next PROPERTY if $property_meta->class_name ne $class;
        my $property_name = $property_meta->property_name;
        next PROPERTY if $property_name eq $target_name_pl_ub;
        my $new_value = $self->$property_name;
        next PROPERTY if not defined $new_value;
        $self->status_message("Property: $property_name");
        $self->status_message('To: '.$self->display_name_for_value($new_value));
        $properties_requested_to_update++;
        OBJECT: for my $obj ( @objects ) {
            if ( grep { $property_name eq $_ } @$only_if_null 
                    and not $property_meta->is_add_remove
                    and defined( my $value = $obj->$property_name) ) {
                my $obj_name = $self->display_name_for_value($obj);
                my $value_name = $self->display_name_for_value($value);
                $self->error_message("Cannot update $obj_name '$property_name' because it already has a value: $value_name");
                next OBJECT;
            }
            my $rv = eval{ $obj->$property_name( $self->$property_name ); };
            if ( defined $rv ) {
                $self->status_message('Successfully updated: '.$self->display_name_for_value($obj));
                $success++; 
            } 
            else { 
                $self->error_message('Failed to update: '.$self->display_name_for_value($obj));
            }
        }
    }

    my $attempted = @objects * $properties_requested_to_update;
    $self->status_message(
        sprintf(
            "Update complete.\nAttempted: %s\nFailed: %s\nSuccess: %s\n",
            $attempted,
            $attempted - $success,
            $success,
        )
    );

    return ( $success ? 1 : 0 ); 
}

1;

