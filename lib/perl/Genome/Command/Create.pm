package Genome::Command::Create;

use strict;
use warnings;

use Genome;
      
require Lingua::EN::Inflect;

class Genome::Command::Create {
    is => 'Command::V2',
    is_abstract => 1,
    doc => 'CRUD create command class.',
};

sub _target_class { Carp::confess('Please use CRUD or implement _target_class in '.$_[0]->class); }
sub _target_name { Carp::confess('Please use CRUD or implement _target_name in '.$_[0]->class); }
sub _target_name_a { return Lingua::EN::Inflect::A($_[0]->_target_name); }
sub _before { return; };

sub sub_command_sort_position { .1 };

sub help_brief { return $_[0]->_target_name_a; }

sub help_detail {
    return "This command creates ".$_[0]->_target_name_a.'.';
}

sub execute {
    my $self = shift;

    $self->status_message('Create '.$self->_target_name.'...');

    $self->_before();

    my $class = $self->class;
    my @properties = grep { $_->class_name eq $class } $self->__meta__->property_metas;
    my (%attrs, $display_string_for_attrs);
    for my $property ( @properties ) {
        my $property_name = $property->property_name;
        my @values = $self->$property_name;
        next if not defined $values[0];
        if ( $property->is_many ) {
            $attrs{$property_name} = \@values;
        }
        else {
            $attrs{$property_name} = $values[0];
        }
        $display_string_for_attrs .= $property_name.' => ',$self->_display_name_for_value($attrs{$property_name});
    }

    $self->status_message($display_string_for_attrs);
    my $target_class = $self->_target_class;
    my $obj = eval { $target_class->create(%attrs) };
    if ( not $obj ) {
        $self->error_message('Could not create ' . $target_class . ": $@");
        return;
    }
    $self->status_message('ID is: '.$obj->id);
    $self->status_message('Create complete. Committing...');

    return 1;
}

sub _display_string_for_attrs {
    my ($self, %attrs ) = @_;

    my $string;
    for my $name ( sort keys %attrs ) {
        $string .= $name.' => ';
        if ( not ref $attrs{$name} or @{$attrs{$name}} <= 10 ) {
            $string .= $self->_display_name_for_value($attrs{$name});
        }
        else {
            $string .= @{$attrs{$name}};
        }
        $string .= "\n";
    }

    return $string;
}

1;

