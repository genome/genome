package Genome::Config::AnalysisProject::SubjectMapping::Input;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::SubjectMapping::Input {
    is => [ "Genome::Utility::ObjectWithTimestamps", "Genome::Utility::ObjectWithCreatedBy" ],
    table_name => 'config.subject_mapping_input',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        subject_mapping => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping',
            id_by => 'subject_mapping_id',
            constraint_name => 'subject_mapping_input_subject_mapping_id_fkey',
        },
        key => { is => 'Text' },
        value => { is => 'Text' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

1;
