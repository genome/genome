package Genome::Model::Event::Output;

#:eclark 11/17/2009 Code review.

# This output tracking class should become part of either SoftwareResult or the Workflow system.

use strict;
use warnings;

use Genome;
class Genome::Model::Event::Output {
    type_name => 'genome model event output',
    table_name => 'model.event_output',
    id_by => [
        event           => { is => 'Genome::Model::Event', id_by => 'event_id', constraint_name => 'GMEO_GME_FK' },
        name            => { is => 'VARCHAR2', len => 100, column_name => 'PARAM_NAME' },
        value           => { is => 'VARCHAR2', len => 1000, column_name => 'PARAM_VALUE' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
