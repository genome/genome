package Genome::Site::TGI::Synchronize::Classes::ProjectInstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::ProjectInstrumentData {
    table_name => <<SQL
    (
        select part.project_id project_id, cast(part.part_id as integer) part_id
        from genome_project_part part
        where part.label = 'instrument_data'
    ) project_part
SQL
,
    id_by => [
        project_id => { is => 'Integer', },
        part_id => { is => 'Integer', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

