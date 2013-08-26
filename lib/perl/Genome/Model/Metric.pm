package Genome::Model::Metric;

use strict;
use warnings;

use Genome;

class Genome::Model::Metric {
    type_name => 'genome model metric',
    table_name => 'model.build_metric',
    id_by => [
        build => { is => 'Genome::Model::Build', id_by => 'build_id', constraint_name => 'GMM_BI_FK' },
        value => { is => 'VARCHAR2', len => 1000, column_name => 'METRIC_VALUE' },
        name  => { is => 'VARCHAR2', len => 100, column_name => 'METRIC_NAME' },
    ],
    has => [
        model          => { is => 'Genome::Model', via => 'build' },
        model_id       => { via => 'build' },
        model_name     => { via => 'model', to => 'name' },
        data_directory => { via => 'build' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;

