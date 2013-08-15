package Genome::Site::TGI::Synchronize::Classes::SetupProjectSample; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Synchronize::Classes::SetupProjectSample {
    table_name => <<SQL
    (
        --Setup Research Projects Sample
        select distinct p.setup_project_id project_id, woi.dna_id sample_id
        from setup_project p
        join gsc.setup_work_order wo on wo.project_id = p.setup_project_id
        join gsc.work_order_item woi on woi.setup_wo_id = wo.setup_wo_id
        where p.setup_project_id > 2570000 and p.project_type = 'setup project research' and woi.dna_id is not null
         union
        --Setup Work Order Sample
        select distinct woi.setup_wo_id project_id, woi.dna_id sample_id
        from work_order_item woi
        where woi.setup_wo_id > 2570000 and woi.dna_id is not null
    ) work_order_sequence_product
SQL
    ,
    id_by => [
        project_id => { is => 'Integer', },
        sample_id => { is => 'Integer', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::Oltp',
};

1;

