package Genome::Model::Command::UnarchiveInstrumentData;

use strict;
use warnings;

use Genome;

use feature qw(say);

class Genome::Model::Command::UnarchiveInstrumentData {
    is => 'Command::V2',
    has_input => [
        models => {
            is => 'Genome::Model',
            doc => 'models to unarchive',
            is_many => 1,
            shell_args_position => 1,
        },
        volume => {
            is => 'Genome::Disk::Volume',
            doc => 'A location to which scratch allocation information can be stored during unarchiving',
            is_optional => 1,
        },
    ],
    doc => 'a helper that produces commands that can be run to unarchive data',
};

sub execute {
    my $self = shift;

    for my $model ($self->models) {
        my $unarchive_cmd = Genome::InstrumentData::Command::Unarchive->create(
            instrument_data => [$model->instrument_data],
            analysis_project => $model->analysis_project,
            volume => $self->volume,
        );
        unless ($unarchive_cmd) {
            $self->fatal_message('Failed to create unarchive instrument data command!');
        }
        unless ($unarchive_cmd->execute) {
            $self->fatal_message('Failed to execute unarchive instrument data command!');
        }
    }

    return 1;
}


1;
