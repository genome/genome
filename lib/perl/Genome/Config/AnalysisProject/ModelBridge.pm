package Genome::Config::AnalysisProject::ModelBridge;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::ModelBridge {
    is => [ "Genome::Utility::ObjectWithTimestamps", "Genome::Utility::ObjectWithCreatedBy" ],
    table_name => 'config.analysis_project_model_bridge',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id',
            constraint_name => 'analysis_project_model_bridge_analysis_project_id_fkey',
        },
        model => {
            is => 'Genome::Model',
            id_by => 'model_id',
            constraint_name => 'analysis_project_model_bridge_model_id_fkey',
        },
        config_profile_item => {
            is => 'Genome::Config::Profile::Item',
            id_by => 'profile_item_id',
            is_optional => 1,
            constraint_name => 'analysis_project_model_bridge_profile_item_id_fkey',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

__PACKAGE__->add_observer(aspect => 'create', callback => \&_is_created);

sub _is_created {
    my $self = shift;
    Genome::Timeline::Event::AnalysisProject->model_created(
        $self->model->id,
        $self->analysis_project,
    );
}

1;
