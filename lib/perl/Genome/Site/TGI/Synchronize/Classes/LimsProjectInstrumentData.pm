package Genome::Site::TGI::Synchronize::Classes::LimsProjectInstrumentData; 

use strict;
use warnings;

use Genome;

# EX: project id 2744972
class Genome::Site::TGI::Synchronize::Classes::LimsProjectInstrumentData {
    table_name => <<SQL
    (
		--Administration Projects Instrument Data
		select distinct ap.project_id project_id, woisp.seq_id instrument_data_id
		from administration_project ap
		join project_work_order pwo on pwo.project_id = ap.project_id
		join work_order_item woi on woi.setup_wo_id = pwo.setup_wo_id
		join gsc.woi_sequence_product\@dw woisp on woisp.woi_id = woi.woi_id
		where ap.project_id > 2570000
		 and ap.status != 'abandoned'
		union
		--Setup Work Order Instrument Data
		select distinct woi.setup_wo_id project_id, woisp.seq_id instrument_data_id
		from setup s
		join work_order_item woi on woi.setup_wo_id = s.setup_id
		join gsc.woi_sequence_product\@dw woisp on woisp.woi_id = woi.woi_id
		where s.setup_id > 2570000
		 and s.setup_status != 'abandoned'
    ) lims_project_instrument_data
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

