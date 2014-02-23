package Genome::Config::AnalysisProject::SubjectMapping;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping {
    is           => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    id_generator => '-uuid',
    data_source  => 'Genome::DataSource::GMSchema',
    table_name   => 'config.subject_mapping',
    id_by        => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        analysis_project => {
            is    => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id'
        },
        subject_bridges => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping::Subject',
            reverse_as => 'subject_mapping',
            is_many => 1,
        },
        subjects => {
            via => 'subject_bridges',
            to => 'subject',
        },
        inputs => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping::Input',
            reverse_as => 'subject_mapping',
            is_many => 1,
        }
    ]
};

1;
