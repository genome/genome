package Genome::Model::Event::Metric;

#:eclark 11/17/2009 Code review.

# This metric tracking class should become part of either SoftwareResult or the Workflow system.

use strict;
use warnings;

use Genome;
class Genome::Model::Event::Metric {
    table_name => 'model.event_metric',
    type_name => 'genome model event metric',
    id_by => [
        event => {
            is => 'Genome::Model::Event',
            id_by => 'event_id',
            constraint_name => 'GMEM_GME_FK',
        },
        name => { is => 'Text', len => 100, column_name => 'METRIC_NAME' },
    ],
    has => [
        value => {
            is => 'Text',
            len => 1000,
            column_name => 'METRIC_VALUE',
            is_optional => 1,
        },
        model => { is => 'Genome::Model', via => 'event' },
        event_type => { via => 'event' },
        model_name => { via => 'model', to => 'name' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

