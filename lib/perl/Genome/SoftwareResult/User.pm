package Genome::SoftwareResult::User;

use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::User {
    table_name => 'result.user',
    id_by => [
        id => { is => 'Text', len => 32 },
    ],
    has => [
        software_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'software_result_id',
            constraint_name => 'SRU_SR_FK',
        },
        user_class => {
            is => 'UR::Object::Type',
            id_by => 'user_class_name',
        },
        user_id => { is => 'Text', len => 256 },
        user => {
            is => 'UR::Object',
            id_by => 'user_id',
            id_class_by => 'user_class_name',
        },
        active => {
            is => 'Boolean',
            len => 1,
            default_value => 1,
            doc => 'Results actively being used should not be deleted',
        },
        label => { is => 'Text' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    doc => 'links a software result to other entities which depend on it',
};

1;

