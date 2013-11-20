package Genome::Site::TGI::Synchronize::Classes::ProjectInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::ProjectInstrumentData {
    table_name => <<SQL
    (
        select part.project_id project_id, part.part_id entity_id
        from subject.project_part part
        where part.label = 'instrument_data'
    ) project_part
SQL
,
    id_by => [
        project_id => { is => 'Text', },
        entity_id => { is => 'Text', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

