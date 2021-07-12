package Genome::Model::Metric;

use strict;
use warnings;

use Genome;

class Genome::Model::Metric {
    table_name => 'model.build_metric',
    type_name => 'genome model metric',
    id_by => [
        build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
            constraint_name => 'GMM_BI_FK',
        },
        value => {
            is => 'Text',
            len => 1000,
            column_name => 'metric_value',
        },
        name => { is => 'Text', len => 100, column_name => 'metric_name' },
    ],
    has => [
        model => { is => 'Genome::Model', via => 'build' },
        model_id => { via => 'build' },
        model_name => { via => 'model', to => 'name' },
        data_directory => { via => 'build' },
        metric_name => { via => '__self__', to => 'name' },
        metric_value => { via => '__self__', to => 'value' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

