package Genome::FeatureList::Command::List;

use strict;
use warnings;

use Genome;

class Genome::FeatureList::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::FeatureList'
        },
        show => { default_value => 'id,name,source,format,content_type,reference'},
    ],
    doc => 'list feature-lists',
};

1;
