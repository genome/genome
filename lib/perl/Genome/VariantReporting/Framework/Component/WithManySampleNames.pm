package Genome::VariantReporting::Framework::Component::WithManySampleNames;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::Component::WithManySampleNames {
    is => ['Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    has => [
        sample_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
        },
    ],
};

1;
