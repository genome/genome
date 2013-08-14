package Genome::Site::TGI::EntityAlias; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::EntityAlias {
    table_name => 'ENTITY_ALIAS',
    id_by => [
        entity_id => { 
            is => 'Text', 
            column_name => 'ENTITY_ID', 
        },
        entity_type_name => { 
            is => 'Text', 
            column_name => 'ENTITY_TYPE_NAME', 
        },
        alias => {
            is => 'Text', 
            column_name => 'ALIAS', 
        },
        alias_source => {
            is => 'Text', 
            column_name => 'ALIAS_SOURCE', 
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

1;

