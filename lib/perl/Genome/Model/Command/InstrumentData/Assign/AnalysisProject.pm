package Genome::Model::Command::InstrumentData::Assign::AnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData::Assign::AnalysisProject {
    is => 'Genome::Model::Command::InstrumentData::Assign::Base',
    has_input => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'analysis project of instrument data to assign',
        },
    ],
};

sub _resolve_instrument_data {
    my $self = shift;
    my $model = $self->model;
    unless($model->can('compatible_instrument_data')) {
        die $self->error_message('Models of type %s do not currently support finding compatible data.');
    }

    my @compatible_instrument_data = $self->model->compatible_instrument_data;
    my $analysis_project_id = $self->analysis_project->id;

    return grep { $_->analysis_project_bridges(analysis_project_id => $analysis_project_id) } @compatible_instrument_data;
}

1;
