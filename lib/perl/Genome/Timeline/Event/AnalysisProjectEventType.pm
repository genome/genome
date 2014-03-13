package Genome::Timeline::Event::AnalysisProjectEventType;

use strict;
use warnings;

use Genome;

class Genome::Timeline::Event::AnalysisProjectEventType {
    is => 'UR::Object',
    table_name => 'timeline.analysis_project_event_type',
    data_source => 'Genome::DataSource::GMSchema',
    id_by => [
        id => {
            is => 'Text',
            valid_values => [
                'model_created',
                'instrument_data_assigned',
                'status_changed',
                'cle_changed',
                'config_added',
            ]
        }
    ],
    has => [
        events => {
            is_many => 1,
            is => 'Genome::Timeline::Event::AnalysisProject',
            reverse_as => 'type',
        }
    ],
};

1;
