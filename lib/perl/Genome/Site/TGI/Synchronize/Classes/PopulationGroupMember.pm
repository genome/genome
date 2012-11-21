package Genome::Site::TGI::Synchronize::Classes::PopulationGroupMember; 

use strict;
use warnings;

=pod
MEMBER_ID NUMBER   (10)                    {null} NOT NULL ok
PG_ID     NUMBER   (10)                    {null} NOT NULL ok [population_group_id]
properties 2, copied 2, updated 0
=cut

class Genome::Site::TGI::Synchronize::Classes::PopulationGroupMember {
    is => 'UR::Object',
    table_name => 'GSC.POPULATION_GROUP_MEMBER',
    id_by => [
        member_id => { is => 'Number', },
        population_group_id => { is => 'Number', },
    ],
    data_source => 'Genome::DataSource::GMSchema',
};

sub properties_to_copy {# 2
    return (qw/
        member_id
        population_group_id
        /);
}

sub properties_to_keep_updated {# 0
    return;
}

1;

