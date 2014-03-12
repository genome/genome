package Genome::Config::AnalysisProject::ModelBridge;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::ModelBridge {
    is => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'config.analysis_project_model_bridge',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id',
        },
        model => {
            is => 'Genome::Model',
            id_by => 'model_id',
        },
    ],
};

1;
