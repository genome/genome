package Genome::MiscAttribute;

use strict;
use warnings;

use Genome;
class Genome::MiscAttribute {
    table_name => 'subject.misc_attribute',
    id_by => [
        entity_id => { is => 'Text', len => 1000 },
        entity_class_name => { is => 'Text', len => 255 },
        property_name => { is => 'Text', len => 255 },
    ],
    has => [
        value => { is => 'Text', len => 4000, is_optional => 1 },
        entity => {
            calculate_from => [ 'entity_id', 'entity_class_name' ],
            calculate => q( return $entity_class_name->get($entity_id) ),
        },
        _instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'entity_id',
        },
        _model => { is => 'Genome::Model', id_by => 'entity_id' },
        _build => { is => 'Genome::Model::Build', id_by => 'entity_id' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
