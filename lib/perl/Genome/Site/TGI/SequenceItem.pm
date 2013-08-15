package Genome::Site::TGI::SequenceItem; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::SequenceItem {
    table_name => 'SEQUENCE_Item',
    id_by => [
        seq_id => { 
            is => 'Number', 
            column_name => 'seq_id', 
        },
    ],
    has => [
        sequence_item_name => {
            is => 'Text',
            column_name => 'sequence_item_name',
        },
        sequence_item_type => {
            is => 'Text',
            column_name => 'sequence_item_type',
        },
        seq_length => {
            is => 'Text',
            column_name => 'seq_length',
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

1;

