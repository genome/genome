package Genome::InstrumentData::Reimport;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Reimport{ 
};

sub attribute_label_for_reimported_from { 'reimported_from' }

sub attribute_labels_to_ignore_when_reimporting {
    (qw/ bam_path genotype_file genotype_file_name ignored import_date import_format user_name /);
}

sub reimported_from {
    my ($self, $instrument_data) = @_;

    Carp::confess('No instrument data to get which instrument data it was reimported from.') if not $instrument_data;

    my $attribute_label_for_reimported_from = $self->attribute_label_for_reimported_from;
    my $attribute = $instrument_data->attributes(attribute_label => $attribute_label_for_reimported_from);
    return if not $attribute;

    return Genome::InstrumentData->get($attribute->attribute_value);
}

1;

