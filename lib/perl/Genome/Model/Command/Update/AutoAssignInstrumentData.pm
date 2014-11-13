package Genome::Model::Command::Update::AutoAssignInstrumentData;

class Genome::Model::Command::Update::AutoAssignInstrumentData {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'Models for which auto_assign_inst_data will be set to the value provided.',
        },
        value => {
            is => 'Text',
            valid_values => ['0', '1'],
            doc => 'Enable or disable the auto_assign_inst_data flag.',
        },
    ],
};

sub help_detail {
    return 'Set auto_assign_inst_data to the value for the models.'
}


sub execute {
    my $self = shift;

    my @models = $self->models;
    my $value = $self->value;

    for my $model (@models) {
        $model->auto_assign_inst_data($value);
    }

    return 1;
}

1;
