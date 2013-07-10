package Genome::Site::TGI::Sample::Attribute;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::Sample::Attribute {
    table_name => 'SAMPLE_ATTRIBUTE',
    id_by => [
        sample_id        => { column_name => 'ORGANISM_SAMPLE_ID' },
        nomenclature    => { default_value => 'WUGC', doc => 'the since an attribute name can have different meanings in different contexts, this supplies context' },
        name            => { column_name => 'ATTRIBUTE_LABEL' },
    ],
    has_optional => [
        value           => { column_name => 'ATTRIBUTE_VALUE' },
        sample          => { is => 'Genome::Site::TGI::Sample', id_by => 'sample_id' },
        sample_name     => { via => 'sample', to => 'name' },
        patient         => { is => 'Genome::Site::TGI::Individual', via => 'sample', to => 'source' },
        common_name     => { via => 'patient', doc => 'patient common name' },
    ],
    data_source => 'Genome::DataSource::Dwrac',
    doc => 'attributes of samples'
};

1;

__END__

