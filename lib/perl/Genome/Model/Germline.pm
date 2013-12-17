package Genome::Model::Germline;

use strict;
use warnings;
use Genome;

class Genome::Model::Germline {
    is => 'Genome::ModelDeprecated',
    has_param => [
        regions_file => {
            type => 'String',
            doc => "Regions File that defines ROI",
        },
    ],
    has => [
        source_model => {
            is => 'Genome::Model::ReferenceAlignment',
            id_by => 'source_model_id',
        },
        source_model_id => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'source_id', value_class_name => 'Genome::Model::ReferenceAlignment' ],
            is_many => 0,
            is_mutable => 1,
        },
        server_dispatch => {
            is_constant => 1,
            is_class_wide => 1,
            value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
            doc => 'lsf queue to submit the launcher or \'inline\''
        },
        job_dispatch => {
            is_constant => 1,
            is_class_wide => 1,
            value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
            doc => 'lsf queue to submit jobs or \'inline\' to run them in the launcher'
        },
    ],
};

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

1;
