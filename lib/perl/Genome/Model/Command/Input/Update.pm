package Genome::Model::Command::Input::Update;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Input::Update {
    is => 'Genome::Model::Command::Input::Base',
    english_name => 'genome model input command add',
    doc => 'Update inputs to a model',
    has => [
        name => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'The name of the input property to update.',
        },
        value => {
            is => 'Text',
            is_optional => 1,
            shell_args_position => 3,
            doc => 'The value to set for the input. To set a property to NULL, do not give a value. A required property cannot be set to NULL.',
        },
    ],
};

sub help_detail {
    return <<HELP;
    This command will update an input from a model. The input must be an singular property, meaning there must be only one value (not 'is_many'). If the property has more one value, use the 'add' or 'remove' commands to modify the input.
HELP
}

sub execute {
    my $self = shift;

    $self->status_message('Update input');

    for my $model ($self->models) {

        my $name = $self->name;
        unless ( $name ) {
            $self->error_message('No name given to update model input.');
            return;
        }

        my @properties = $model->real_input_properties;
        return if not @properties;
        my ($property) = grep { $name eq $_->{name} } @properties;
        if( not $property ) {
            #if no property found, try falling back on the name of the input itself
            ($property) = grep { $name eq $_->{input_name} } @properties;
            $name = $property->{name} if $property;
        }

        if ( not $property ) {
            $self->error_message("Failed to find input property for $name. Here are the valid input names:");
            $self->show_inputs_and_values_for_model($model);
            return;
        }

        if ( $property->{is_many} ) {
            $self->error_message("Cannot use this update command for a property ($name) with many values. Use the add or remove command.");
            return;
        }

        $self->status_message('Model: '.$model->__display_name__);
        $self->status_message('Property: '.$name);

        my $value = $self->value;
        if ( not defined $value ) {
            if ( not $property->{is_optional} ) {
                $self->error_message("Cannot update $name to NULL because it is a required model input!");
                return;
            }
            # FIXME - do we want to ask here?
            #$self->status_message("Really update '$name' to NULL/undefined? (y/[n])");
            #my $response = <STDIN>;
            #return if $response !~ /y/i;
            $self->status_message('Value: NULL');

            my $input_name = $property->{input_name};

            if( defined $model->$name() ) { #only remove if there's a value
                my $rv = $model->remove_input(name => $input_name);
                unless($rv eq 1) {
                    $self->error_message('Failed to update!');
                    return;
                }
            }
        }
        else {
            $self->status_message('Value filter: '.$value);
            $value = $self->values_from_property_for_filter($property, $value);
            return if not defined $value;
            my $display_name = $self->display_name_for_value($value);
            $self->status_message('Value display name: '.$display_name);
            my $rv = $model->$name($value);
            my $rv_display_name = $self->display_name_for_value($rv);
            if ( not defined $rv or $display_name ne $rv_display_name ) {
                $self->error_message('Failed to update!');
                return;
            }
        }
    }
    $self->status_message('Update successful!');

    return 1;
}

1;

