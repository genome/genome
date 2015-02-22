package Genome::VariantReporting::Framework::Component::WithManyLibraryNames;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::Component::WithManyLibraryNames {
    is => ['Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    has => [
        library_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
            doc => 'List of library names to be used in the report',
        },
        library_name_labels => {
            is => 'HASH',
            is_translated => 1,
            is_optional => 1,
            default => {},
            doc => 'Hash of library_name to label',
        }
    ],
};

sub create_library_specific_field_name {
    my ($self, $field, $library_name) = @_;

    return join("_", $self->get_library_label($library_name), $field);
}

sub get_library_label {
    my ($self, $library_name) = @_;

    return $self->library_name_labels->{$library_name} || $library_name;
}

sub create_library_specific_field_names {
    my ($self, $fields) = @_;

    my @field_names;
    for my $library_name ($self->library_names) {
        for my $field (@$fields) {
            push @field_names, $self->create_library_specific_field_name($field, $library_name);
        }
    }
    return @field_names;
}
1;

