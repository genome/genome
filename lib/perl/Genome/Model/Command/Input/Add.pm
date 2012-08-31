package Genome::Model::Command::Input::Add;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Input::Add {
    is => 'Genome::Model::Command::Input::Base',
    english_name => 'genome model input command add',
    doc => 'Add inputs to a model',
    has => [
        name => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'The name of the input property to add.',
        },
        'values' => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 3,
            doc => 'The values of the inputs. Id or resolved via filter string. Ex: library_id=$ID'
        },
    ],
};

sub help_detail {
    return <<EOS;
    This command will add inputs from a model. The input must be an 'is_many' property, meaning there must be more than one input allowed (eg: instrument_data). If the property is singular, use the 'update' command.  Use the plural name of the property.
EOS
}

sub execute {
    my $self = shift;

    $self->status_message('Add inputs');

    for my $model ($self->models) {
        if ( not $model ) {
            $self->error_message('No model to update input');
            return;
        }

        my $name = $self->name;
        if ( not $name ) {
            $self->error_message('No name given to update model input.');
            return;
        }

        my @properties = $model->real_input_properties;
        return if not @properties;
        my ($property) = grep { $name eq $_->{name} } @properties;
        if ( not $property ) {
            $self->error_message("Failed to find input property for $name. Here are the valid input names:");
            $self->show_inputs_and_values_for_model($model);
            return;
        }

        if ( not $property->{is_many} ) {
            $self->error_message("Cannot use this add command for a property ($name) with many values. Use the update command.");
            return;
        }

        my @values = $self->values;
        if ( not @values ) {
            $self->error_message('No values to add');
            return;
        }

        $self->status_message('Model: '.$model->__display_name__);
        $self->status_message('Property: '.$name);
        $self->status_message('Value filters: '.join(' ', @values));

        # get values for filters
        my %values = $self->unique_values_from_property_for_filter($property, @values);
        return if not %values;

        # get existing values
        my %existing_values;
        for my $value ( $model->$name ) {
            my $id = $self->display_name_for_value($value);
            $existing_values{$id} = $value;
        }

        # add
        my $add_method = $property->{add_method};
        for my $value_id ( keys %values ) {
            if ( exists $existing_values{$value_id} ) {
                $self->status_message('Value exists on model, skipping: '.$value_id);
                next;
            }
            $self->status_message('Add: '.$value_id);
            $model->$add_method($values{$value_id});
            $existing_values{$value_id} = $values{$value_id};
        }

        # verify
        for my $value ( $model->$name ) { 
            my $id = $self->display_name_for_value($value);
            delete $values{$id};
        }
        if ( %values ) {
            $self->error_message('Failed to add these values: '.join(' ', keys %values));
            return;
        }
    }
    $self->status_message('Add successful!');

    return 1; 
}

1;

