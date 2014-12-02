package Genome::VariantReporting::Framework::Test::WithTranslationsInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::WithTranslationsInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    has_translated => {
        translated1 => {},
        translated2 => {
            is_many => 1,
        },
        translated3 => {
            is_many => 1,
        },
    },
};

sub name {
    return '__translated_test__';
}

sub requires_annotations {
    return qw(__translated_test__);
}

sub field_descriptions {
    return (
        info => 'Description of info',
    );
}

sub _interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele}->{info} = pp($entry->info_for_allele($variant_allele));
    }

    return %return_values;
}

1;
