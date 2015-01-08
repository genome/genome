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
            doc => 'List of sample names to be used in the report',
        },
        sample_name_labels => {
            is => 'HASH',
            is_translated => 1,
            is_optional => 1,
            default => {},
            doc => 'Hash of sample_name to label',
        }
    ],
};

sub create_sample_specific_field_name {
    my ($self, $field, $sample_name) = @_;

    return join("_", $self->get_sample_label($sample_name), $field);
}

sub get_sample_label {
    my ($self, $sample_name) = @_;


    return $self->sample_name_labels->{$sample_name} || $sample_name;
}

sub create_sample_specific_field_names {
    my ($self, $fields) = @_;

    my @field_names;
    for my $sample_name ($self->sample_names) {
        for my $field (@$fields) {
            push @field_names, $self->create_sample_specific_field_name($field, $sample_name);
        }
    }
    return @field_names;
}

1;
