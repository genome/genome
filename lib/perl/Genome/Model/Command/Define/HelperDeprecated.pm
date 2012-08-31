package Genome::Model::Command::Define::HelperDeprecated;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::HelperDeprecated {
    is => 'Genome::Model::Command::Define::Helper',
    is_abstract => 1,
    has => [
        processing_profile_name => {
            is => 'Text',
            via => 'processing_profile',
            to => 'name',
        },
    ],
    has_optional => [
        subject_name => {
            is => 'Text',
            via => 'subject',
            to => 'name',
            doc => 'The name of the subject of the model',
        },
    ],
};

1;
