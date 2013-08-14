package Genome::SoftwareResult::Metric;

use strict;
use warnings;

use Genome;
class Genome::SoftwareResult::Metric {
    table_name => 'result.metric',
    type_name => 'software result metric',
    id_by => [
        metric_name => {
            is => 'VARCHAR2',
            len => 1000,
        },
        software_result_id => {
            is => 'Text',
            len => 32,
        },
    ],
    has => [
        metric_value => {
            is => 'VARCHAR2',
            len => 1000,
        },
        software_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'software_result_id',
            constraint_name => 'SRM_SR_FK',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
