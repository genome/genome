package Genome::Site::TGI::Synchronize::Classes::LimsProjectInstrumentData; 

use strict;
use warnings;

use Genome;

# EX: project id 2744972
class Genome::Site::TGI::Synchronize::Classes::LimsProjectInstrumentData {
    table_name => <<SQL
    (
		--Setup Research Projects Seq Product
		select distinct s.setup_id project_id, woisp.seq_id instrument_data_id
		from administration_project ap
		join project_work_order pwo on pwo.project_id = ap.project_id
		join work_order_item woi on woi.setup_wo_id = pwo.setup_wo_id
		join gsc.woi_sequence_product\@dw woisp on woisp.woi_id = woi.woi_id
		join setup s on s.setup_name = ap.project_name
        join setup_project sp on sp.setup_project_id = s.setup_id
		where s.setup_id > 2570000 and sp.project_type = 'setup project research'
		union
		--Setup Work Order Seq Product
		select distinct woi.setup_wo_id project_id, woisp.seq_id instrument_data_id
		from gsc.woi_sequence_product\@dw woisp
		join work_order_item woi on woi.woi_id = woisp.woi_id
		join setup_work_order wo on wo.setup_wo_id = woi.setup_wo_id
		where woi.setup_wo_id > 2570000
      ) work_order_sequence_product
SQL
    ,
    id_by => [
        project_id => { is => 'Text', },
        instrument_data_id => { is => 'Text', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::Oltp',
};

1;

