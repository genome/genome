# review gsanders 
# I know not much about this subject but I do not understand why this is so similar to variantreviewlistmember.pm...
# Are these redundant or am I missing the relationship?

package Genome::VRListMember;

use strict;
use warnings;

use Genome;
class Genome::VRListMember {
    type_name => 'genome vr list member',
    table_name => 'GENOME_VR_LIST_MEMBER',
    id_by => [
        #list_id   => { is => 'NUMBER', len => 10 },
        list => { is => 'Genome::VRList', id_by => 'list_id', },
        #member_id => { is => 'NUMBER', len => 10 },
        member => { is => 'Genome::VariantReviewDetail', id_by => 'member_id', },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
