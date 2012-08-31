package Genome::Model::Command::Update::ApipeCronStatus;

class Genome::Model::Command::Update::ApipeCronStatus {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'Models for which apipe_cron_status will be set. Resolved by Genome::Command::Base.',
        },
        value => {
            is => 'Text',
            doc => 'The value to use for the APipe cron status.',
        },
    ],
};

sub _is_hidden_in_doc { return !Genome::Sys->current_user_is_admin };

sub help_detail {
    return 'Sets the APipe Cron Status note for the models to the value provided.'
}


sub execute {
    my $self = shift;

    unless (Genome::Sys->current_user_is_admin) {
        print "You must be an admin to execute this command!\n";
        return;
    }

    my @models = $self->models;
    my $value = $self->value;

    for my $model (@models) {
        $model->set_apipe_cron_status($value);
    }

    return 1;
}

1;
