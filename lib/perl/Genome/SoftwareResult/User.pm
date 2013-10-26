package Genome::SoftwareResult::User;

use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::User {
    table_name => 'result.user',
    id_generator => '-uuid',
    id_by => [
        id => {
            is => 'Text',
            len => 32,
        },
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
        user_id => {
            is => 'VARCHAR2',
            len => 256,
        },
        user => {
            is => 'UR::Object',
            id_by => 'user_id',
            id_class_by => 'user_class_name',
        },
        active => {
            is => 'BOOLEAN',
            doc => 'Results actively being used should not be deleted',
            default => 1,
        },
        label => {  },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'links a software result to other entities which depend on it',
};

1;

