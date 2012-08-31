package Genome::Model::Tools::DetectVariants2::Result::Filter;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Filter {
    is => ['Genome::Model::Tools::DetectVariants2::Result::DetectionBase'],
    has_param => [
        filter_name => {
            is => 'Text',
            doc => 'The name of the filter to use',
        },
        filter_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Additional parameters to pass to the filter',
        },
        filter_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the filter to use',
        },
        previous_filter_strategy => {
            is => 'Text',
            is_optional => 1,
            doc => 'the dispatcher string corresponding to the previously run filters for the data being filtered',
        },
    ],
};

#Most filter-specific logic is in Detector.pm

1;
