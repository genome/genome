package Genome::Model::Command::SplitOnLibrary;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::SplitOnLibrary {
    is => ['Command::V2'],
    has => [
        model => {
            is => 'Genome::Model',
            doc => 'model to split',
            shell_args_position => 1,
        },
    ],
    has_optional_transient_output => [
        new_models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'the newly created models',
        },
    ],
    doc => 'command to make per-library models based on a per-sample model',
};

sub help_detail {
    return <<HELP
This command copies the provided model once for each library of the instrument data assigned to it, and assigns only that instrument data for that library.
HELP
}

sub execute {
    my $self = shift;

    my $model = $self->model;

    unless($model->can('instrument_data')) {
        die $self->error_message("Model cannot have instrument_data assigned");
    }

    my @instrument_data = $model->instrument_data;
    unless(@instrument_data) {
        die $self->error_message("Model has no instrument_data assigned");
    }

    my @libraries = Genome::Library->get([map $_->library_id, @instrument_data]);
    unless(scalar(@libraries) > 0) {
        die $self->error_message("Error retrieving libraries for instrument data.");
    }

    if(scalar(@libraries) == 1) {
        die $self->error_message("Instrument data from one library... nothing to split.");
    }

    my @new_models;
    for my $library (@libraries) {
        push @new_models, $self->_make_per_library_model($model, $library);
    }

    $self->status_message("New models: %s", join(", ", map $_->id, @new_models));
    $self->new_models(\@new_models);

    return 1;
}

sub _make_per_library_model {
    my $self = shift;
    my $template = shift;
    my $library = shift;

    my @instrument_data_to_assign = grep { $_->library_id eq $library->id } $template->instrument_data;
    my $new_model_name = sprintf("%s-%s", $template->name, $library->name);
    my $model = $template->copy(instrument_data => \@instrument_data_to_assign, name => $new_model_name);
    $model->auto_assign_inst_data(0);

    return $model;
}

1;
