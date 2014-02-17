package Genome::Config::AnalysisProject::SubjectMapping::Subject;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Subject {
    is           => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    id_generator => '-uuid',
    data_source  => 'Genome::DataSource::GMSchema',
    table_name   => 'config.subject_mapping_subject',
    id_by        => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        subject => {
            is    => 'Genome::Subject',
            id_by => 'subject_id'
        },
        label => {
            is => 'Text',
        }
    ]
};

1;
