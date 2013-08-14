package Genome::Site::TGI::Synchronize::Classes::ProjectSample;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::ProjectSample {
    table_name => <<SQL
    (
        select part.project_id project_id, cast(part.part_id as integer) part_id
        from genome_project_part part
        where part.label = 'sample'
    ) project_part
SQL
,
    id_by => [
        project_id => { is => 'Text', },
        part_id => { is => 'Text', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

