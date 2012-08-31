package Genome::SampleSource;

use strict;
use warnings;

use Genome;

class Genome::SampleSource {
    is => ['Genome::Subject','Genome::Searchable'],
    is_abstract => 1,
    has => [
        taxon_id => {
            is => 'Number',
            via => 'attributes',
            to => 'attribute_value',
            where => [ attribute_label => 'taxon_id' ],
            is_mutable => 1,
            doc => 'Taxon ID for this source',
        },
        taxon => {
            is => 'Genome::Taxon',
            id_by => 'taxon_id',
        },
        species_name => {
            via => 'taxon'
        },
    ],
};

1;

