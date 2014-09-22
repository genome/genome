package Genome::Config::AnalysisProject::SubjectMapping::Subject;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Subject {
    is => [ "Genome::Utility::ObjectWithTimestamps", "Genome::Utility::ObjectWithCreatedBy" ],
    table_name => 'config.subject_mapping_subject',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        subject_mapping => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping',
            id_by => 'subject_mapping_id',
            constraint_name => 'subject_mapping_subject_subject_mapping_id_fkey',
        },
        subject => {
            is => 'Genome::Subject',
            id_by => 'subject_id',
            constraint_name => 'subject_mapping_subject_subject_id_fkey',
        },
        label => { is => 'Text' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

1;
