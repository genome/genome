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

sub create_sample_specific_field_name {
    my ($self, $field, $sample_name) = @_;

    return join("_", $sample_name, $field);
}

sub create_sample_specific_field_names {
    my ($self, $fields, $sample_names) = @_;

    my @field_names;
    for my $sample_name (@$sample_names) {
        for my $field (@$fields) {
            push @field_names, $self->create_sample_specific_field_name($field, $sample_name);
        }
    }
    return @field_names;
}

1;
