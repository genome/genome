package Genome::Model::Command::InstrumentData::Unassign;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData::Unassign {
    is => 'Genome::Command::Base',
    has => [
        model => {
            is => 'Genome::Model', 
            shell_args_position => 1,
            doc => 'Model resolved via expression.'
        },
       instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            is_optional => 1,
            doc => 'Instrument data to remove resolved via expression.',
        },
        all => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'Remove all instrument data from the model.',
        },
    ],
    doc => "remove instrument data from a model",
};

sub execute {
    my $self = shift;

    $self->status_message('Remove instrument data...');

    my $model = $self->model;
    $self->status_message('Model: '. $model->__display_name__);

    my @instrument_data = $self->instrument_data;
    if ( $self->all ) {
        if ( @instrument_data > 0 ) {
            $self->error_message('Cannot remove all and specific instrument data!');
            return;
        }
        @instrument_data = $model->instrument_data;
    }

    for my $instrument_data ( @instrument_data ) {
        if ( $model->remove_instrument_data($instrument_data) ) {
            $self->status_message($instrument_data->id.' removed');
        }
        else {
            $self->status_message($instrument_data->id.' not found on model');
        }
    }

    $self->status_message('Remove instrument data...OK');

    return 1;
}

1;

