package Genome::VariantReviewListFilter;

use strict;
use warnings;

use Genome;
class Genome::VariantReviewListFilter {
    type_name => 'variant review list filter',
    table_name => 'VARIANT_REVIEW_LIST_FILTER',
    id_by => [
        filter_name => { is => 'VARCHAR2', len => 256 },
        list_id     => { is => 'NUMBER', len => 10 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
