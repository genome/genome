package Genome::Command::WithColor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Command::WithColor {
    is => ['Command::V2', 'Genome::Utility::ColorMixin'],
    has => [
        color => {
            is => 'Boolean',
            is_param => 1,
            is_optional => 1,
            default_value => 1,
            doc => 'Use colors in display.'
        },
    ],
};

1;
