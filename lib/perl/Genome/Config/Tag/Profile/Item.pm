package Genome::Config::Tag::Profile::Item;

use strict;
use warnings;

use Genome;

class Genome::Config::Tag::Profile::Item {
    is => [ "Genome::Utility::ObjectWithTimestamps", "Genome::Utility::ObjectWithCreatedBy" ],
    table_name => 'config.tag_profile_item',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        tag => {
            is => 'Genome::Config::Tag',
            id_by => 'tag_id',
            constraint_name => 'tag_profile_item_tag_id_fkey',
        },
        tag_name => { is => 'Text', via => 'tag', to => 'name' },
        profile_item => {
            is => 'Genome::Config::Profile::Item',
            id_by => 'profile_item_id',
            constraint_name => 'tag_profile_item_profile_item_id_fkey',
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            via => 'profile_item',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

1;
