package Genome::Config::Tag;

use strict;
use warnings;

use Genome;

class Genome::Config::Tag {
    is => [ "Genome::Utility::ObjectWithTimestamps", "Genome::Utility::ObjectWithCreatedBy" ],
    table_name => 'config.tag',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        name => {
            is => 'Text',
            doc => 'The canonical "name" of the tag',
        },
        description => {
            is => 'Text',
            is_optional => 1,
            doc => 'Free-form explanatory information about the tag',
        },
        profile_item_bridges => {
            is => 'Genome::Config::Tag::Profile::Item',
            reverse_as => 'tag',
            is_optional => 1,
            is_many => 1,
        },
        profile_items => {
            is => 'Genome::Config::Profile::Item',
            via => 'profile_item_bridges',
            to => 'profile_item',
            is_optional => 1,
            is_many => 1,
            is_mutable => 1,
        },
        subject_mapping_bridges => {
            is => 'Genome::Config::Tag::AnalysisProject::SubjectMapping',
            reverse_as => 'tag',
            is_optional => 1,
            is_many => 1,
        },
        subject_mappings => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping',
            via => 'subject_mapping_bridges',
            to => 'subject_mapping',
            is_optional => 1,
            is_many => 1,
            is_mutable => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

1;
