package Genome::Config::Tag::AnalysisProject::SubjectMapping;

use strict;
use warnings;

use Genome;

class Genome::Config::Tag::AnalysisProject::SubjectMapping {
    is => [ "Genome::Utility::ObjectWithTimestamps", "Genome::Utility::ObjectWithCreatedBy" ],
    table_name => 'config.tag_subject_mapping',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        tag => {
            is => 'Genome::Config::Tag',
            id_by => 'tag_id',
            constraint_name => 'tag_subject_mapping_tag_id_fkey',
        },
        tag_name => { is => 'Text', via => 'tag', to => 'name' },
        subject_mapping => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping',
            id_by => 'subject_mapping_id',
            constraint_name => 'tag_subject_mapping_subject_mapping_id_fkey',
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            via => 'subject_mapping',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

1;
