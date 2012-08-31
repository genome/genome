package Genome::SoftwareResult::User;

use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::User {
    table_name => 'SOFTWARE_RESULT_USER',
    id_by => [
        id => { is => 'NUMBER', len => 10 },
    ],
    has => [
        software_result     => { is => 'Genome::SoftwareResult', id_by => 'software_result_id' },
        user_class          => { is => 'UR::Object::Type', id_by => 'user_class_name' },
        user_id             => { is => 'VARCHAR2', len => 256 },
        user                => { is => 'UR::Object', id_class_by => 'user_class_name', id_by => 'user_id' },
        label               => { }
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc         => 'links a software result to other entities which depend on it',
};

1;

