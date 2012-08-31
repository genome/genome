package Genome::SoftwareResult::Metric;

use strict;
use warnings;

use Genome;
class Genome::SoftwareResult::Metric {
    type_name => 'software result metric',
    table_name => 'SOFTWARE_RESULT_METRIC',
    id_by => [
        metric_name        => { is => 'VARCHAR2', len => 100 },
        software_result_id => { is => 'NUMBER', len => 20 },
    ],
    has => [
        metric_value                    => { is => 'VARCHAR2', len => 1000 },
        software_result                 => { is => 'Genome::SoftwareResult', id_by => 'software_result_id', constraint_name => 'SRM_SR_FK' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
