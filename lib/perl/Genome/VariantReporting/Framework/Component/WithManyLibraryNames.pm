package Genome::VariantReporting::Framework::Component::WithManyLibraryNames;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::Component::WithManyLibraryNames {
    is => ['Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    has_transient_optional => [
        library_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
            doc => 'List of library names to be used in the report',
        },
    ],
};

sub create_library_specific_field_name {
    my ($self, $field, $library_name) = @_;

    return join("_", $library_name, $field);
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

