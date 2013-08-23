package Genome::Site::TGI::Synchronize::Classes::SetupProjectSequenceProduct; 

use strict;
use warnings;

use Genome;

# EX: project id 2744972
class Genome::Site::TGI::Synchronize::Classes::SetupProjectSequenceProduct {
    table_name => <<SQL
    (
        --Setup Research Projects Seq Product
        select distinct p.setup_project_id project_id, woisp.seq_id
        from setup_project p
        join gsc.setup_work_order wo on wo.project_id = p.setup_project_id
        join gsc.work_order_item woi on woi.setup_wo_id = wo.setup_wo_id
        join gsc.woi_sequence_product\@dw woisp on woisp.woi_id = woi.woi_id
        where p.setup_project_id > 2570000 and p.project_type = 'setup project research'
        union
        --Setup Work Order Seq Product
        select distinct woi.setup_wo_id project_id, woisp.seq_id seq_id
        from gsc.woi_sequence_product\@dw woisp
        join work_order_item woi on woi.woi_id = woisp.woi_id
        join setup_work_order wo on wo.setup_wo_Id = woi.setup_wo_id
        where woi.setup_wo_id > 2570000
    ) work_order_sequence_product
SQL
    ,
    id_by => [
        project_id => { is => 'Integer', },
        seq_id => { is => 'Integer', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::Oltp',
};

1;

