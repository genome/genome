package Genome::Config::AnalysisProject::SubjectMapping::Input;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Input {
    is           => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    id_generator => '-uuid',
    data_source  => 'Genome::DataSource::GMSchema',
    table_name   => 'config.subject_mapping_input',
    id_by        => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        subject_mapping => {
            is    => 'Genome::Config::AnalysisProject::SubjectMapping',
            id_by => 'subject_mapping_id'
        },
        key => {
            is => 'Text',
        },
        value => {
            is => 'Text',
        }
    ]
};

1;
