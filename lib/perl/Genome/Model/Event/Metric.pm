package Genome::Model::Event::Metric;

#:eclark 11/17/2009 Code review.

# This metric tracking class should become part of either SoftwareResult or the Workflow system.

use strict;
use warnings;

use Genome;
class Genome::Model::Event::Metric {
    type_name => 'genome model event metric',
    table_name => 'model.event_metric',
    id_by => [
        event   => { is => 'Genome::Model::Event', id_by => 'event_id', constraint_name => 'GMEM_GME_FK' },
        name    => { is => 'VARCHAR2', len => 100, column_name => 'METRIC_NAME' },
    ],
    has => [
        value   => { is => 'VARCHAR2', len => 1000, is_optional => 1, column_name => 'METRIC_VALUE' },
        model   => { is => 'Genome::Model', via => 'event' },
        event_type  => { via => 'event', to => 'event_type' },
        model_name  => { via => 'model', to => 'name' },
        
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

