package Genome::VariantReporting::Framework::Test::WithTranslationsFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::WithTranslationsFilter {
    is => 'Genome::VariantReporting::Framework::Component::Filter',
    has_translated => {
        translated2 => {},
    },
};

sub name {
    return '__translated_test_filter__';
}

sub requires_annotations {
    return qw(__translated_test__);
}

1;
