package Genome::Model::MetagenomicCompositionShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicCompositionShotgun {
    is => 'Genome::ModelDeprecated',
    has => [
    map({
            $_ => {
                via => 'processing_profile',
            }
        } Genome::ProcessingProfile::MetagenomicCompositionShotgun->params_for_class
    ),
    ],
};
# This model type has 'from_model_links'. When this model was deleted, these 'sub-models' were also deleted.

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

sub build_subclass_name {
    return 'metagenomic-composition-shotgun';
}

1;

