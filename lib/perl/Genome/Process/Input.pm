package Genome::Process::Input;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Process::Input {
    table_name => 'process.input',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    id_by => [
        id => {is => 'Text'},
    ],
    has => [
        process => {
            is => 'Genome::Process',
            id_by => 'process_id',
        },
        process_id => {
            is => 'Text',
        },
        name => {
            is => 'Text',
            doc => 'The name of the input property on the process',
        },
        value_class_name => {
            is => 'Text',
            doc => 'The class_name of the object this input refers to.',
        },
        value_id => {
            is => 'Text',
            doc => 'The id of the object this input refers to.',
        },
        value => {
            is_optional => 1,
            is => 'UR::Object',
            id_by => 'value_id',
            id_class_by => 'value_class_name',
        },
        array_index => {
            is => 'Number',
            default => -1,
        }
    ],
};


1;
