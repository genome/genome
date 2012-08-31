package Genome::Model::Command::Input::Remove;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Input::Remove {
    is => 'Genome::Model::Command::Input::Base',
    english_name => 'genome model input command remove',
    doc => 'Remove inputs to a model.',
    has => [
        name => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'The name of the input property to remove.',
        },
        'values' => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 3,
            doc => 'The values of the inputs. Id or resolved via filter string. Ex: library_id=$ID'
        },
    ],
    has_optional => [
        abandon_builds => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Abandon (not remove) builds that have these inputs.',
        },
    ],
};

sub help_detail {
    return <<EOS;
    This command will remove inputs from a model. The input must be an 'is_many' property, meaning there must be more than one input allowed (eg: instrument_data). If the property is singular, use the 'update' command.  Use the plural name of the property.
    
    Optionally, builds associated with these inputs may be abandoned.  Note that the builds are not removed, only abandoned.
EOS
}
sub execute {
    my $self = shift;

    $self->status_message('Remove inputs');

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
            $self->error_message("Cannot use this remove command for a property ($name) with many values. Use the update command.");
            return;
        }

        my @values = $self->values;
        if ( not @values ) {
            $self->error_message('No values to remove');
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
            my $id = ( $value->can('id') ? $value->id : $value );
            $existing_values{$id} = $value;
        }

        # remove
        my $transaction = UR::Context::Transaction->begin;
        my $remove_method = $property->{remove_method};
        for my $value_id ( keys %values ) {
            if ( not defined $existing_values{$value_id} ) {
                $self->status_message('Value does not exist on model, skipping: '.$value_id);
                next;
            }
            $self->status_message('Remove: '.$value_id);
            $model->$remove_method($values{$value_id});
            $existing_values{$value_id} = $values{$value_id};
        }

        # verify
        my @existing_values = $model->$name;
        if ( not $property->{is_optional} and not @existing_values ) {
            $transaction->rollback;
            $self->error_message('Cannot remove all values for required property!');
            return;
        }
        my @failed_to_remove;
        for my $value ( @existing_values ) { 
            my $id = ( $value->can('id') ? $value->id : $value );
            push @failed_to_remove, $id if defined $values{$id};
        }
        if ( @failed_to_remove ) {
            $self->error_message('Failed to remove these values: '.join(' ', @failed_to_remove));
            return;
        }
        $transaction->commit;

        # abandon
        if ( $self->abandon_builds ) {
            $self->_abandon_builds( keys %values );
        }
    }
    $self->status_message('Remove successful!');

    return 1; 
}

sub _abandon_builds {
    my ($self, @value_ids) = @_;

    my @builds = Genome::Model::Build->get(model_id => [map($_->id, $self->models)]);
    if ( not @builds ) {
        $self->status_message('No builds to abandon. Model does not have any');
        return 1;
    }

    my $name = $self->name;
    my @builds_to_abandon;
    BUILD: for my $build ( @builds ) {
        my @values = $build->$name;
        next if not @values;
        VALUE: for my $value ( @values ) {
            my $id = $self->display_name_for_value($value);
            next VALUE if not grep { $id eq $_ } @value_ids;
            push @builds_to_abandon, $build;
            next BUILD;
        }
    }

    if ( not @builds_to_abandon ) {
        $self->status_message('No builds to abandon. There are no builds with the removed inputs');
        return 1;
    }

    for my $build ( @builds ) {
        my $rv = eval{ $build->abandon; }; # this can die
        if ( $rv ) {
            $self->status_message('Abandon build: '.$build->__display_name__);
        }
        else {
            $self->status_message('Failed to abandon build: '.$build->__display_name__);
        }
    }

    return 1;
}

1;

