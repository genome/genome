package Genome::Model::Command::Update::BuildRequested;

class Genome::Model::Command::Update::BuildRequested {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'Models for which build_requested will be set to the value provided. Resolved by Genome::Command::Base.',
        },
        value => {
            is => 'Text',
            valid_values => ['-1', '0', '1'],
            doc => 'Enable or disable the build_requested flag.',
        },
    ],
};

sub _is_hidden_in_docs { return !Genome::Sys->current_user_is_admin };

sub help_detail {
    return 'Set build_requested to the value for the models.'
}


sub execute {
    my $self = shift;

    unless (Genome::Sys->current_user_is_admin) {
        print "You must be an admin to use this command.\n";
        return;
    }

    my @models = $self->models;
    my $value = $self->value;

    for my $model (@models) {
        $model->build_requested($value);
    }

    return 1;
}

1;
