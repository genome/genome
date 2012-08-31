# review gsanders 
# I know not much about this subject but I do not understand why this is so similar to variantreviewlist.pm...
# Are these redundant or am I missing the relationship?

package Genome::VRList;

use strict;
use warnings;

use Genome;
class Genome::VRList {
    type_name => 'genome vr list',
    table_name => 'GENOME_VR_LIST',
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
    has_many => [
        members => { 
            is => 'Genome::VRListMember',
            reverse_id_by => 'list', 
        },
        details => { 
            is => 'Genome::VariantReviewListDetail',
            via => 'members',
            to => 'member',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

1;
