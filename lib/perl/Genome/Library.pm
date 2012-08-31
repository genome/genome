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
        fragment_size_range => { 
            is => 'Text', 
            column_name => 'LIBRARY_INSERT_SIZE',
            doc => 'intended size range of fragments from library construction' 
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
        protocol_name => { 
            is_transient => 1, 
            is => 'Text', 
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
        # what is this??? -ss
        protocol_name           => { is_transient => 1, is => 'Text', },
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

