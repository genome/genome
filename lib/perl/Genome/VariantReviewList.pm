# review gsanders 
# I know not much about this subject but I do not understand why this is so similar to VRList.pm...
# Are these redundant or am I missing the relationship?

package Genome::VariantReviewList;

use strict;
use warnings;

use Genome;
class Genome::VariantReviewList {
    type_name => 'variant review list',
    table_name => 'VARIANT_REVIEW_LIST',
    id_by => [
        id => { is => 'NUMBER', len => 10 },
    ],
    has => [
        author         => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        name           => { is => 'VARCHAR2', len => 256, is_optional => 1 },
        position_lists => { is => 'VARCHAR2', len => 1024, is_optional => 1 },
        rt_ticket      => { is => 'NUMBER', len => 10, is_optional => 1 },
        source_maps    => { is => 'VARCHAR2', len => 1024, is_optional => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
