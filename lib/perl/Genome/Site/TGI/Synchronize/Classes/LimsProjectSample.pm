package Genome::Site::TGI::Synchronize::Classes::LimsProjectSample; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::LimsProjectSample {
    table_name => <<SQL
    (
        --Administration/Setup Research Project Sample
        select distinct s.setup_id project_id, wois.organism_sample_id sample_id
        from administration_project ap
		join project_work_order pwo on pwo.project_id = ap.project_id
		join work_order_item woi on woi.setup_wo_id = pwo.setup_wo_id
		join woi_sample wois on wois.woi_id = woi.woi_id
		join setup s on s.setup_name = ap.project_name
        join setup_project sp on sp.setup_project_id = s.setup_id
		where s.setup_id > 2570000 and sp.project_type = 'setup project research'
		union
		--Setup Work Order Sample
		select distinct woi.setup_wo_id project_id, wois.organism_sample_id sample_id
		from work_order_item woi
		join woi_sample wois on wois.woi_id = woi.woi_id
		where woi.setup_wo_id > 2570000
    ) lims_project_sample
SQL
    ,
    id_by => [
        project_id => { is => 'Text', },
        sample_id => { is => 'Text', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::Oltp',
};

1;

