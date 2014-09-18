package Genome::InstrumentData::Reimport;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Reimport{ 
};

sub attribute_labels_to_ignore_when_reimporting {
    (qw/ bam_path genotype_file genotype_file_name
        ignored import_date import_format is_paired_end
        original_data_path original_data_path_md5
        read_length read_count user_name /);
}

sub attributes_for_reimport_from_instrument_data {
    my ($self, $instrument_data) = @_;

    die 'No instrument data given!' if not $instrument_data;

    my $reimported_attribute_label = $self->attribute_label_for_reimported_from;
    my %reimport = ( 
        $reimported_attribute_label => $instrument_data->id,
        library_name => $instrument_data->library->name,
    );

    my $source_file = eval{ $instrument_data->bam_path; };
    if ( not $source_file ) {
        $source_file = eval{ $instrument_data->archive_path; };
    }
    if ( not $source_file or not -s $source_file ) {
        $self->error_message("Source file for instrument data (%s) does not exist!", $instrument_data->id);
        return;
    }
    $reimport{source_files} = $source_file;

    for my $optional_property_name (qw/ run_name subset_name /) {
        next if not $instrument_data->$optional_property_name;
        $reimport{$optional_property_name} = $instrument_data->$optional_property_name;
    }

    ATTRIBUTE: for my $attribute ( $instrument_data->attributes ) {
        for my $attribute_label_to_ignore ( $self->attribute_labels_to_ignore_when_reimporting ) {
            next ATTRIBUTE if $attribute->attribute_label eq $attribute_label_to_ignore;
        }
        $reimport{ $attribute->attribute_label } = $attribute->attribute_value;
    }

    return \%reimport;
}

sub headers_for_reimport_attributes {
    my ($self, @reimports) = @_;

    Carp::confess('No reimports to get headers!') if not @reimports;

    my %headers = map { $_ => 1 } map { keys %$_ } @reimports;
    for (qw/ library_name source_files /) { delete $headers{$_}; }
    my @headers = sort keys %headers;
    unshift @headers, (qw/ library_name source_files /);

    return @headers;
}

sub attribute_label_for_reimported_from { 'reimported_from' }

sub reimported_from {
    my ($self, $instrument_data) = @_;

    Carp::confess('No instrument data to get which instrument data it was reimported from.') if not $instrument_data;

    my $attribute_label_for_reimported_from = $self->attribute_label_for_reimported_from;
    my $attribute = $instrument_data->attributes(attribute_label => $attribute_label_for_reimported_from);
    return if not $attribute;

    return Genome::InstrumentData->get($attribute->attribute_value);
}

1;

