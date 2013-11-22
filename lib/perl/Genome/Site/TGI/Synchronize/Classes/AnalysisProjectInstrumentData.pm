package Genome::Site::TGI::Synchronize::Classes::AnalysisProjectInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::AnalysisProjectInstrumentData {
    table_name => <<SQL
    (
        select instrument_data_id, analysis_project_id
        from config.instrument_data_analysis_project_bridge
    )
    ap_instdata
SQL
    ,
    data_source => 'Genome::DataSource::GMSchema',
    id_by => [
        instrument_data_id => { is => 'Text', },
        analysis_project_id => { is => 'Text', },
    ],
};

1;

