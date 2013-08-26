package Genome::Task::Params;

use strict;
use warnings;

use Command::Dispatch::Shell;
use Genome;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use JSON::XS;

class Genome::Task::Params {
    table_name => 'web.task_params',
    id_by => [
        genome_task_id => {
            is => 'Text',
            len => 255,
        },
    ],
    has => [
        task => {
            is => 'Genome::Task',
            id_by => 'genome_task_id',
            constraint_name => 'GTP_FK',
        },
        params => {
            is=>'Text',
            doc => 'JSON encoded param hash',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'params for scheduled tasks',
};

1;
