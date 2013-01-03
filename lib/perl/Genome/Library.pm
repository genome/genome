package Genome::Library;

use strict;
use warnings;
use Genome;

class Genome::Library {
    is => ['Genome::Notable','Genome::Searchable'],
    id_by => [
        library_id => {
            is => 'Number',
        },
    ],
    has => [
        name => {
            is => 'Text',
            column_name => 'FULL_NAME',
            doc => 'Name of the library. Usually has the sample name and an extension.',
        },
        sample_id => {
            is => 'Text',
        },
        sample => {
            is => 'Genome::Sample',
            id_by => 'sample_id',
            doc => 'Sample that this library came from.',
        },
        sample_name => {
            is => 'Text',
            via => 'sample',
            to => 'name'
        },
    ],
    has_optional => [
        original_insert_size => {
            is => 'Text',
            column_name => 'ORIGINAL_INSERT_SIZE',
            doc => 'The original insert size of the fragments. This may be a number or a range.'
        },
        library_insert_size => {
            is => 'Text',
            column_name => 'LIBRARY_INSERT_SIZE',
            doc => 'The relative insert size of fragments. This may be a number or a range.'
        },
        fragment_size_range => {
            is => 'Text',
            column_name => 'LIBRARY_INSERT_SIZE',
        },
        protocol => {
            is => 'Text',
            column_name => 'PROTOCOL',
            doc => 'Protocol used to generate the library.'
        },
        taxon_id => {
            is => 'Number',
            via => 'sample',
        },
        taxon => {
            is => 'Genome::Taxon',
            via => 'sample',
        },
        species_name => {
            is => 'Text',
            via => 'taxon',
        },
        sample_source => {
            via => 'sample',
            to => 'source',
            doc => 'Source of the sample',
        },
        sample_source_name => {
            via => 'sample_source',
            to => 'name',
            doc => 'Name of the sample\'s source'
        },
        sample_source_id => {
            via => 'sample_source',
            to => 'id',
            doc => 'ID of the sample\'s source'
        },
        models => {
            is => 'Genome::Model',
            via => 'sample',
            to => 'models',
            is_many => 1,
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            reverse_as => 'library',
            is_many => 1,
        },
    ],
    table_name => 'FRAGMENT_LIBRARY',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub __display_name__ {
    return $_[0]->name.' ('.$_[0]->id.')';
}

sub delete {
    my $self = shift;

    $self->status_message("Deleting library " . $self->__display_name__);

    my @instrument_data = Genome::InstrumentData->get(
        library_id => $self->id,
    );
    for my $instrument_data (@instrument_data) {
        $instrument_data->delete;
    }

    return $self->SUPER::delete(@_);
}

1;

