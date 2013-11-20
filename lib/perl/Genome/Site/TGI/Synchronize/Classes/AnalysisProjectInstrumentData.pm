package Genome::Site::TGI::Synchronize::Classes::AnalysisProjectInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::AnalysisProjectInstrumentData {
    is => 'UR::Object',
    table_name => 'config.instrument_data_analysis_project_bridge',
    data_source => 'Genome::DataSource::GMSchema',
    id_by => [
        instrument_data_id => { is => 'Text', },
        analysis_project_id => { is => 'Text', },
    ],
};

1;

