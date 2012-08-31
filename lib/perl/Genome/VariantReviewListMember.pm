# review gsanders 
# I know not much about this subject but I do not understand why this is so similar to VRListMember.pm...
# Are these redundant or am I missing the relationship?

package Genome::VariantReviewListMember;

use strict;
use warnings;

use Genome;
class Genome::VariantReviewListMember {
    type_name => 'variant review list member',
    table_name => 'VARIANT_REVIEW_LIST_MEMBER',
    id_by => [
        list_id   => { is => 'NUMBER', len => 10 },
        member_id => { is => 'NUMBER', len => 10 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
