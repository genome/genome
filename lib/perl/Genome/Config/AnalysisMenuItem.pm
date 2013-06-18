package Genome::Config::AnalysisMenuItem;

use strict;
use warnings;

class Genome::Config::AnalysisMenuItem {
    is => 'Genome::Utility::ObjectWithTimestamps',
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'GENOME_ANALYSIS_MENU_ITEM',
    id_by => [
        id => {
            is => 'Text',
        },
    ],
    has => [
        name => {
            is => 'Text',
        },
        configuration_set_id => {
            is => 'Text',
        },
        path => {
            is => 'Text',
            via => 'configuration_set',
            to => 'path',
        },
        configuration_set => {
            is => 'Genome::Config::Set',
            id_by => 'configuration_set_id',
        },
        created_at => {
            is => 'Timestamp',
        },
        updated_at => {
            is => 'Timestamp',
        },
    ],
};

1;
