package Genome::Site::TGI::Synchronize::Classes::LimsProjectSample; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsProjectSample {
    is => 'Genome::Site::TGI::Synchronize::Classes::LimsProjectPartBase',
    table_name => <<SQL
    (
        --Administration Project Sample
		select distinct ap.project_id project_id, wois.organism_sample_id entity_id
		from administration_project ap
		join project_work_order pwo on pwo.project_id = ap.project_id
		join work_order_item woi on woi.setup_wo_id = pwo.setup_wo_id
		join woi_sample wois on wois.woi_id = woi.woi_id
		where ap.project_id > 2570000
         and ap.status != 'abandoned'
		union
		--Setup Work Order Sample
		select distinct woi.setup_wo_id project_id, wois.organism_sample_id entity_id
		from setup s
		join work_order_item woi on woi.setup_wo_id = s.setup_id
		join woi_sample wois on wois.woi_id = woi.woi_id
		where woi.setup_wo_id > 2570000
         and s.setup_status != 'abandoned'
    ) lims_project_sample
SQL
    ,
    id_by => [
        project_id => { is => 'Text', },
        entity_id => { is => 'Text', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::Oltp',
};

1;

