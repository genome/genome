package Genome::Site::TGI::WorkOrderItemSequenceProduct; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::WorkOrderItemSequenceProduct {
    table_name => 'WOI_SEQUENCE_PRODUCT',
    id_by => [
        work_order_item => {
            is => 'Genome::WorkOrderItem',
            id_by => 'woi_id',
        },
        sequence_item => {
            is => 'Genome::Site::TGI::SequenceItem',
            id_by => 'seq_id',
        },
    ],
    has_optional => [ 
        instrument_data => { # optional b/c seq item may be a read or entire region
            is => 'Genome::InstrumentData',
            id_by => 'seq_id',
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

1;

