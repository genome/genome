package Genome::Site::TGI::Synchronize::Classes::SampleAttribute;

use strict;
use warnings;


class Genome::Site::TGI::Synchronize::Classes::SampleAttribute {
    table_name => 'SAMPLE_ATTRIBUTE',
    id_by => [
        organism_sample_id => {
            is => 'Number',
        },
        attribute_label => {
            is => 'Text',
        },
        attribute_value => {
            is => 'Text',
        },
        nomenclature => {
            is => 'Text',
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

1;

