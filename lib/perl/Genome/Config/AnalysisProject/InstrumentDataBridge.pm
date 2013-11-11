package Genome::Config::AnalysisProject::InstrumentDataBridge;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::InstrumentDataBridge {
    is => ['Genome::Utility::ObjectWithTimestamps'],
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'config.instrument_data_analysis_project_bridge',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        analysis_project => {
            is    => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id'
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'instrument_data_id',
        },
        status => {
            is => 'Text',
            valid_values => ['new','failed','processed'],
        },
        fail_count => {
            is => 'Integer',
            default_value => 0,
        },
        reason => {
            is => 'Text'
        },
    ],
};


1;
