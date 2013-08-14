package Genome::Site::TGI::WorkOrderItem; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::WorkOrderItem {
    table_name => <<SQL
    (
        select woi_id id, setup_wo_id setup_id from work_order_item
    ) work_order_item
SQL
    ,
    id_by => [
        id => {
            is => 'Text'
        },
    ],
    has=> [ 
        setup_id => {
            is => 'Text'
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

1;
