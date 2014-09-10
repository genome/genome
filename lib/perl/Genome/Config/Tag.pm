package Genome::Config::Tag;

use strict;
use warnings;

use Genome;

class Genome::Config::Tag {
    is => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'config.tag',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        name => {
            is => 'Text',
            doc => 'The canonical "name" of the tag',
        },
        description => {
            is_optional => 1,
            doc => 'Free-form explanatory information about the tag',
            is => 'Text',
        },
        profile_item_bridges => {
            is => 'Genome::Config::Tag::Profile::Item',
            is_many => 1,
            is_optional => 1,
            reverse_as => 'tag',
        },
        profile_items => {
            is => 'Genome::Config::Profile::Item',
            via => 'profile_item_bridges',
            to => 'profile_item',
            is_many => 1,
            is_optional => 1,
        },
    ],
};

1;
