package Genome::Site::TGI::Synchronize::Classes::SampleAttribute;

use strict;
use warnings;


class Genome::Site::TGI::Synchronize::Classes::SampleAttribute {
    table_name => 'GSC.SAMPLE_ATTRIBUTE',
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
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

